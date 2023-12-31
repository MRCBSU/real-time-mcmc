suppressMessages(library(lubridate))
suppressMessages(library(tidyverse))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

if(!exists("date.data"))
    date.data <- (today() - days(1)) %>% format("%Y%m%d")

## Where to find the data, if NULL use command line argument
if(!exists("deaths.loc"))
    ## deaths.loc <- paste(date.data, "COVID19 Deaths.csv") ## NULL
	deaths.loc <- paste0("Dataset Modelling " , date.data, ".csv") ## NULL
    #deaths.loc <- file.path(proj.dir, "data/raw/deaths", paste0(date.data, ".csv"))
	#deaths.loc <- paste0("/data/covid-19/data-raw/deaths/", ymd(date.data), ".csv")

## Define an age-grouping
if(!exists("age.agg")){
    age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
    age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+")
    }
nA <- length(age.labs)

if(!exists("regions")){
    regions <- c("London", "Outside_London")
}

## Map our names for columns (LHS) to data column names (RHS)
col.names <- list(
	death_date = "dod",
	finalid = "finalid",
	onset_date = "onsetdate",
	nhs_region = "NHSEnglandRegionName2019_lkup",
	phe_region = "phec_name",
	utla_name = "utla_name"
)

# Given a row in a deaths file, return its region.
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
if(!exists("data.files"))
    data.files <- build.data.filepath("RTM_format/deaths",
                                      "deaths",
                                      date.data,
                                      "_",
                                      regions,
                                      "_",
                                      nA,
                                      "ages.txt")

## Which columns are we interested in?
death.col.args <- list()
death.col.args[[col.names[["death_date"]]]] <- col_character()
death.col.args[[col.names[["onset_date"]]]] <- col_character()
death.col.args[[col.names[["finalid"]]]] <- col_double()
death.col.args[[col.names[["nhs_region"]]]] <- col_character()
death.col.args[[col.names[["phe_region"]]]] <- col_character()
death.col.args[[col.names[["utla_name"]]]] <- col_character()
death.cols <- do.call(cols, death.col.args)	# Calling with a list so use do.call

within.range <- function(dates) {
	return(dates <= today() & dates >= ymd("2020-01-01"))
}

plausible.death.date <- function(x) {
	death.within.range <- within.range(x$Date)
	onset.within.range <- is.na(x$Onset) | within.range(x$Onset)
	after.onset <- is.na(x$Onset) | (x$Onset <= x$Date)
	return(death.within.range & onset.within.range & after.onset)

## Some patients are known to have the month and day swapped
## Fix these here...
fix.dates <- function(df) {
	return(mutate(df,
		orig_date = Date,
		orig_onset = Onset,
		Date = heuristically.swap.day.and.month(Date),
		Onset = heuristically.swap.day.and.month(Onset)
	  ))
}

## Read the file and rename columns
if (is.null(deaths.loc)) {
	input.loc = commandArgs(trailingOnly = TRUE)[1]
} else {
	if (startsWith(deaths.loc, "/")) input.loc <- deaths.loc
	else input.loc = build.data.filepath(subdir = "raw/deaths", deaths.loc)
}
print(paste("Reading from", input.loc))
dth.dat <- read_csv(input.loc,
                    col_types = death.cols) %>%
    rename(!!!col.names) %>%
    mutate(Date = fuzzy_date_parse(death_date),
           Onset = fuzzy_date_parse(onset_date)) %>%
    fix.dates %>%
    mutate(plausible_death_date = plausible.death.date(.) & !is.na(death_date))

if (!all(dth.dat$plausible_death_date)) {
	implausible.dates <- dth.dat %>% filter(!plausible_death_date)
	print("WARNING: The following rows have implausible death dates and have been excluded: ")
	(x.out <- implausible.dates %>% select(c(finalid, Date, Onset))) %>% print(n=1000)
	dth.dat <- dth.dat %>% filter(plausible_death_date)
}

print("WARNING: the following rows have had onset and/or death date changed to become plausible: ")
dth.dat %>%
   filter(orig_date != Date | orig_onset != Onset) %>%
   select(c(finalid, Date, orig_date, Onset, orig_onset)) %>%
   mutate(
   		orig_date = as_date(ifelse(orig_date == Date,
							  NA, orig_date)),
   		orig_onset = as_date(ifelse(orig_onset == Onset,
							  NA, orig_onset))
	) %>%
   print(n=1000)

## The following code was necessary for the first time on 20200324. Hopefully it can be commented out and ignored in future iterations
#dth.dat$dod <- lubridate::as_date(dth.dat$dod, format = "%d/%m/%Y", tz = "GMT")
#dth.dat$dod_NHSE <- lubridate::as_date(dth.dat$dod_NHSE, format = "%m/%d/%Y", tz = "GMT")
#x <- dth.dat$dod
#x[is.na(x)] <- dth.dat$dod_NHSE[is.na(x)]
#dth.dat <- dth.dat %>%
    #mutate(Date = x)
## ## 

latest.date <- ymd(date.data) - reporting.delay
earliest.date <- ymd("2020-02-17")

dth.dat <- dth.dat %>%
    filter(Date <= latest.date) %>%
    filter(Date >= earliest.date) %>%
    filter(age > 0) %>%
    mutate(Region = get.region(.)) %>%
    filter(Region %in% regions) %>%
    mutate(Age.Grp = cut(age, age.agg, age.labs, right = FALSE, ordered_result = T))

if(flg.confirmed)
    dth.dat <- dth.dat %>% filter(death_type == "Confirmed")

rtm.dat <- dth.dat %>%
	group_by(Date, Region, Age.Grp, .drop = FALSE) %>%
	tally %>%
	right_join(		# Add missing rows
		expand_grid(Date = as_date(earliest.date:latest.date), Region = regions, Age.Grp = age.labs),
		by = c("Date", "Region", "Age.Grp")
	) %>%
	replace_na(list(n = 0)) %>%		# 0 if just added
	arrange(Date)

## Write rtm.dat to data file
names(data.files) <- regions
for(reg in regions) {
    region.dat <- pivot_wider(rtm.dat %>%
                              filter(Region == reg),
                              id_cols = 1,
                              names_from = Age.Grp,
                              values_from = n)
    tmpFile <- data.files[reg]
	dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)
    
    print(paste(
        "Writing to",
        tmpFile,
        "(", sum(region.dat[, -1]), "total deaths,", nrow(region.dat), "rows.)"
    ))
    
    region.dat %>%
        write_tsv(
            tmpFile,
            col_names = FALSE
        )
}

## Save the data as processed
write_csv(rtm.dat, file.path(out.dir, "deaths_data.csv"))

## Save a quick plot of the data...
require(ggplot2)
rtm.dat %>%
    group_by(Date, Region) %>%
    summarise(count = sum(n)) -> rtm.dat.plot

gp <- ggplot(rtm.dat.plot, aes(x = Date, y = count, color = Region)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    ggtitle(paste("Daily number of deaths by day of death (on", lubridate::as_date(date.data), ")")) +
    xlab("Date of death") +
    ylab("#Deaths") +
    theme(
        legend.position = "top",
        legend.justification = "left",
        )
plot.filename <- build.data.filepath("RTM_format/deaths", "deaths_plot", date.data, ".pdf")
if (!file.exists(dirname(plot.filename))) dir.create(dirname(plot.filename))
ggsave(plot.filename,
       gp,
       width = 8.15,
       height = 6)
