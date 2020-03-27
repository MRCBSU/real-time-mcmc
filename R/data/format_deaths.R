library(readr)
library(lubridate)
library(dplyr)

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

## YYYYMMDD string, used in filenames and reporting lag
date.data <- "20200325"

## Where to find the data, if NULL use command line argument
linelist.loc <- NULL
# linelist.loc <- paste0(date.data, " - Anonymised Line List.csv")		# relative to data/raw

## Map our names for columns (LHS) to data column names (RHS)
col.names <- list(
	death_date = "dod",
	finalid = "finalid",
	onset_date = "onsetdate"
)

## Inputs that are dependent on the output required.
reporting.delay <- 2
## region.def.str <- "ifelse(nhsregionname == \"London\", \"London\", \"Outside_London\")"
region.def.str <- "\"ENGLAND\""
## death.col.name <- "PATIENT_DEATH_DATE"


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
file.loc <- dirname(thisFile())
proj.dir <- dirname(dirname(file.loc))
dir.data <- file.path(proj.dir, "data")
source(file.path(file.loc, "utils.R"))

## Which columns are we interested in?
death.col.args <- list()
death.col.args[[col.names[["death_date"]]]] <- col_character()
death.col.args[[col.names[["onset_date"]]]] <- col_character()
death.col.args[[col.names[["finalid"]]]] <- col_double()
death.cols <- do.call(cols, death.col.args)	# Calling with a list so use do.call

plausible.death.date <- function(x) {
	within.range <- x$Date <= today() & x$Date >= ymd("2020-01-01")
	after.onset <- is.na(x$Onset) | (x$Onset <= x$Date)
	return(within.range & after.onset)
}

## Some patients are known to have the month and day swapped
## Fix these here...
fix.dates <- function(df) {
	finalid.swapped.onset <- c(2507, 2727, 3159, 3284, 3767, 6969, 35367, 39214)
	finalid.swapped.death <- c(2507)
	df <- df %>% mutate(
		swap_onset = finalid %in% finalid.swapped.onset,
		swap_death = finalid %in% finalid.swapped.death
	)
	df[df$swap_death,]$Date  = swap.day.and.month(df[df$swap_death,]$Date)
	df[df$swap_onset,]$Onset = swap.day.and.month(df[df$swap_onset,]$Onset)
	return(df)
}

## Read the file and rename columns
if (is.null(linelist.loc)) {
	input.loc = commandArgs(trailingOnly = TRUE)[1]
} else {
	input.loc = build.data.filepath(subdir = "raw", linelist.loc)
}
print(paste("Reading from", input.loc))
dth.dat <- read_csv(
		input.loc,
		col_types = death.cols
	) %>%
	rename(!!!col.names) %>%
    mutate(Date = fuzzy_date_parse(death_date),
		   Onset = fuzzy_date_parse(onset_date)) %>%
	fix.dates %>%
	mutate(plausible_death_date = plausible.death.date(.))

if (!all(dth.dat$plausible_death_date)) {
	implausible.dates <- dth.dat %>% filter(!plausible_death_date)
	print("WARNING: The following rows have implausible death dates and have been excluded: ")
	print(implausible.dates %>% select(c(finalid, Date, Onset, swap_death, swap_onset)))
	dth.dat <- dth.dat %>% filter(plausible_death_date)
}

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
all.dates <- as.character(seq(earliest.date, latest.date, by = 1))

dth.dat <- dth.dat %>%
    filter(Date <= latest.date) %>%
    filter(Date >= earliest.date) %>%
    mutate(fDate = factor(Date)) %>%
    mutate(Region = as.factor(eval(parse(text = region.def.str))))
levels(dth.dat$fDate) <- c(levels(dth.dat$fDate), all.dates[!(all.dates %in% levels(dth.dat$fDate))])

rtm.dat <- dth.dat %>%
    group_by(fDate, Region, .drop = FALSE) %>%
    summarise(count = n())

rtm.dat$fDate <- lubridate::as_date(rtm.dat$fDate)

rtm.dat <- arrange(rtm.dat, fDate)

## Write rtm.dat to data file
for(reg in levels(rtm.dat$Region)){
    write.table(filter(rtm.dat, Region == reg) %>%
               select(fDate, count),
            file = build.data.filepath("RTM_format", "deaths", date.data, "_", reg, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
}
