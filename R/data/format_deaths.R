library(assertr)
library(readr)
library(lubridate)
library(dplyr)

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

## Map our names for columns (LHS) to data column names (RHS)
col.names <- list(
	death_date = "dod",
	finalid = "finalid"
)

## YYYYMMDD string, used in filenames and reporting lag
date.data <- "20200325"

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
death.col.args[[col.names[["finalid"]]]] <- col_double()
death.cols <- do.call(cols, death.col.args)	# Calling with a list so use do.call

## Read the file and rename columns
dth.dat <- read_csv(
		build.data.filepath(subdir = "raw", date.data, " COVID19 Deaths.csv"),
		col_types = death.cols
	) %>%
	rename(!!!col.names) %>%
    mutate(Date = fuzzy_date_parse(death_date)) %>%
	## Remove two lines with implausible death dates
	filter(finalid != 9425 & finalid != 30609) %>%
	## Check plausibility (should flag parsing errors)
	verify(Date >= ymd("2020-01-01")) %>%
	verify(Date <= today())

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
