suppressMessages(library(lubridate))
suppressMessages(library(tidyverse))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

## Where to find the data, if NULL use command line argument
deaths.loc <- paste0("/data/covid-19/data-raw/deaths-scotland/", ymd(date.data), ".csv")

## Map our names for columns (LHS) to data column names (RHS)
col.names <- list(
	death_date = "NRS.Date.Death",
	age = "NRS.Age"
)

## Inputs that are dependent on the output required.
reporting.delay <- 3


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
if (!exists("proj.dir")) {
	file.loc <- dirname(thisFile())
	proj.dir <- dirname(dirname(file.loc))
}
dir.data <- file.path(proj.dir, "data")
source(file.path(proj.dir, "R/data/utils.R"))
source(file.path(proj.dir, "config.R"))

## Which columns are we interested in?
death.col.args <- list()
death.col.args[[col.names[["death_date"]]]] <- col_character()
death.col.args[[col.names[["age"]]]] <- col_integer()
death.cols <- do.call(cols_only, death.col.args)	# Calling with a list so use do.call

within.range <- function(dates) {
	return(dates <= today() & dates >= ymd("2020-01-01"))
}

plausible.death.date <- function(x) {
	death.within.range <- within.range(x$Date)
	return(death.within.range)
}

heuristically.swap.day.and.month <- function(dates) {
	swapped.dates <- suppressWarnings(swap.day.and.month(dates))
	should.swap <- (!is.na(swapped.dates) &
					!within.range(dates) &
					within.range(swapped.dates))
	dates[should.swap] <- swapped.dates[should.swap]
	return(dates)
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
input.loc <- deaths.loc
print(paste("Reading from", input.loc))
dth.dat <- read_csv(input.loc,
                    col_types = death.cols) %>%
    rename(!!!col.names) %>%
    mutate(Date = fuzzy_date_parse(death_date)) %>%
    fix.dates %>%
    mutate(plausible_death_date = plausible.death.date(.) & !is.na(death_date),
    	   age_group = cut(age, age.agg, age.labs, right = FALSE, ordered_result = T))

if (!all(dth.dat$plausible_death_date)) {
	implausible.dates <- dth.dat %>% filter(!plausible_death_date)
	print("WARNING: The following rows have implausible death dates and have been excluded: ")
	implausible.dates %>% print(n=1000)
	dth.dat <- dth.dat %>% filter(plausible_death_date)
}

print("WARNING: the following rows have had death dates changed to become plausible: ")
dth.dat %>%
   filter(orig_date != Date) %>%
   select(c(Date, orig_date)) %>%
   mutate(
   		orig_date = as_date(ifelse(orig_date == Date,
							  NA, orig_date))
	) %>%
   print(n=1000)

latest.date <- ymd(date.data) - reporting.delay
earliest.date <- start.date

rtm.dat <- dth.dat %>%
    filter(Date <= latest.date) %>%
    filter(Date >= earliest.date) %>%
	group_by(Date, age_group, .drop = FALSE) %>%
	tally %>%
	right_join(		# Add missing rows
		expand_grid(
			Date = as_date(earliest.date:latest.date),
			age_group = age.labs,
		),
		by = c("Date", "age_group")
	) %>%
	replace_na(list(n = 0))

region.dat <- pivot_wider(
		  rtm.dat,
		  id_cols = 1,
		  names_from = age_group,
		  values_from = n
	)

output.file <- data.files["Scotland"]
print(paste(
	"Writing to",
	output.file,
	"(", sum(rtm.dat$n), "total deaths,", nrow(region.dat), "rows.)"
))

region.dat %>%
	arrange(Date) %>%
	write_tsv(
		output.file,
		col_names = FALSE
	)

## Save a quick plot of the data...
require(ggplot2)
gp <- ggplot(rtm.dat, aes(x = Date, y = n)) +
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
plot.filename <- build.data.filepath("RTM_format/deaths_Scotland", "deaths_plot", date.data, ".pdf")
if (!file.exists(dirname(plot.filename))) dir.create(dirname(plot.filename))
ggsave(plot.filename,
       gp,
       width = 8.15,
       height = 6)
