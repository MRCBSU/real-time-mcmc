suppressMessages(library(lubridate))
suppressMessages(library(readxl))
suppressMessages(library(tidyverse))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

lubridate.date <- ymd(date.data)

# Given a vector, return the first element where predicate is true
# Returns NULL if not true for any
first.where.true <- function(vec, predicate) {
	true.at <- which(predicate(vec))
	if (length(true.at) == 0) return(NULL)
	index.to.use <- min(true.at)
	return(vec[index.to.use])
}	

## Where to find the data, if NULL use command line argument
if(!exists("cases")) {
	possible.cases.locations <- c(
		glue::glue("/data/covid-19/data-raw/dstl/{lubridate.date}/Anonymised Combined Line List {format(lubridate.date, '%Y%m%d')}.xlsx")
	)
	cases.loc <- first.where.true(possible.cases.locations, file.exists)
	if (is.null(cases.loc)) {
		stop(paste('No valid cases data files, tried:', possible.cases.locations))
	}
}


## Map our names for columns (LHS) to data column names (RHS)
possible.col.names <- list(
    case_date = "specimen_date",
    nhs_region = c("NHSER_name", "nhser_name"),
    phe_region = c("PHEC_name", "phec_name"),
    utla_name = c("UTLA_name", "utla_name"),
    age = "age",
    cat = "cat",
    pillar = "pillar"
)
input.col.names <- suppressMessages(names(read_excel(cases.loc, n_max=0)))
is.valid.col.name <- function(name) {name %in% input.col.names}
first.valid.col.name <- function(names) {first.where.true(names, is.valid.col.name)}
col.names <- lapply(possible.col.names, first.valid.col.name)
invalid.col.names <- sapply(col.names, is.null)
if (any(invalid.col.names)) {
	names.invalid.cols <- paste0(names(possible.col.names)[invalid.col.names], collapse = ", ")
	stop(paste("No valid column name for:", names.invalid.cols))
}

# Given a row in a cases file, return its region.
# Various useful functions for this are defined above.
get.region <- nhs.region


####################################################################
## BELOW THIS LINE SHOULD NOT NEED EDITING
####################################################################

## Location of this script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else if (.Platform$GUI == "RStudio" || Sys.getenv("RSTUDIO") == "1") {
                # We're in RStudio
                return(rstudioapi::getSourceEditorContext()$path)
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

## Where are various directories?
if(!exists("file.loc")){
    file.loc <- dirname(thisFile())
    proj.dir <- dirname(dirname(file.loc))
    dir.data <- file.path(proj.dir, "data")
	source(file.path(file.loc, "utils.R"))
	source(file.path(proj.dir, "config.R"))
}

## Which columns are we interested in?
case.cols <- c("numeric", rep("text", 10), "numeric", rep("date", 3), "text")	# Calling with a list so use do.call

within.range <- function(dates) {
	return(dates <= today() & dates >= ymd("2020-01-01"))
}

plausible.case.date <- function(x) {
	return(within.range(x$Date))
}

## Some patients are known to have the month and day swapped
## Fix these here...
fix.dates <- function(df) {
	return(mutate(df,
		orig_date = Date,
		Date = heuristically.swap.day.and.month(Date)
	  ))
}

## Read the file and rename columns
input.loc <- cases.loc
print(paste("Reading from", input.loc))
case.dat <- read_excel(input.loc,
                    col_types = case.cols) %>%
    rename(!!!col.names) %>%
    mutate(Date = fuzzy_date_parse(case_date) %>% na_if(ymd("1990-01-01"))) %>%
    fix.dates %>%
    mutate(plausible_case_date = plausible.case.date(.) & !is.na(case_date))

if (!all(case.dat$plausible_case_date)) {
	implausible.dates <- case.dat %>% filter(!plausible_case_date)
	print("WARNING: The following rows have implausible case dates and have been excluded: ")
	(x.out <- implausible.dates %>% select(c(Date))) %>% print(n=1000)
	case.dat <- case.dat %>% filter(plausible_case_date)
}

print("WARNING: the following rows have had their case date changed to become plausible: ")
case.dat %>%
   filter(orig_date != Date) %>%
   select(Date, orig_date) %>%
   mutate(
   		orig_date = as_date(ifelse(orig_date == Date,
							  NA, orig_date))
	) %>%
   print(n=1000)

latest.date <- ymd(date.data) ## - reporting.delay
earliest.date <- ymd("2020-06-01")

case.dat <- case.dat %>%
    filter(
		Date <= latest.date,
    	Date >= earliest.date,
		cat == "Residential dwelling (including houses, flats, sheltered accommodation)",
		pillar == "Pillar 2"
	) %>%
    mutate(Region = get.region(.)) %>%
    filter(Region %in% regions) %>%
    mutate(Age.Grp = cut(age, age.agg, age.labs, right = FALSE, ordered_result = T))

rtm.dat <- case.dat %>%
	group_by(Date, Region, Age.Grp, .drop = FALSE) %>%
	tally %>%
	right_join(		# Add missing rows
		expand_grid(Date = as_date(start.date:latest.date), Region = regions, Age.Grp = age.labs),
		by = c("Date", "Region", "Age.Grp")
	) %>%
	replace_na(list(n = 0)) %>%		# 0 if just added
	arrange(Date)

## Write rtm.dat to data file
names(cases.files) <- regions
for(reg in regions) {
    region.dat <- pivot_wider(rtm.dat %>%
                              filter(Region == reg),
                              id_cols = 1,
                              names_from = Age.Grp,
                              values_from = n)
    tmpFile <- cases.files[reg]
	dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)
    
    print(paste(
        "Writing to",
        tmpFile,
        "(", sum(region.dat[, -1]), "total cases,", nrow(region.dat), "rows.)"
    ))
    
    region.dat %>%
        write_tsv(
            tmpFile,
            col_names = FALSE
        )
}

## Save the data as processed
save(case.dat, rtm.dat, file = file.path(out.dir, "cases_data.RData"))

## Save a quick plot of the data...
require(ggplot2)
rtm.dat %>%
    group_by(Date, Region) %>%
    summarise(count = sum(n)) %>%
    mutate(ignore = !(Date <= (latest.date - reporting.delay))) -> rtm.dat.plot

gp <- ggplot(rtm.dat.plot, aes(x = Date, y = count, color = Region)) +
    geom_line(aes(linetype = ignore)) +
    geom_point() +
    theme_minimal() +
    ggtitle(paste("Daily number of cases by day of specimen (on", lubridate::as_date(date.data), ")")) +
    xlab("Date of specimen") +
    ylab("#Deaths") +
    theme(
        legend.position = "top",
        legend.justification = "left",
        )
plot.filename <- build.data.filepath("RTM_format/cases", "cases_plot", date.data, "_", reporting.delay, "d.pdf")
if (!file.exists(dirname(plot.filename))) dir.create(dirname(plot.filename))
ggsave(plot.filename,
       gp + guides(linetype=FALSE),
       width = 1.5*8.5,
       height = 1.5*6)

## rtm.dat.plot %>%
##     group_by(Date, ignore) %>%
##     summarise(count = sum(count)) -> rtm.dat.Eng.plot
