library(readr)
library(lubridate)
library(dplyr)

## Get deaths data

## Inputs that should change on a daily basis.
date.data <- "20200324"
dir.data <- "../../../Data/"

## Inputs that are dependent on the output required.
reporting.delay <- 2
## region.def.str <- "ifelse(nhsregionname == \"London\", \"London\", \"Outside_London\")"
region.def.str <- "\"ENGLAND\""
## death.col.name <- "PATIENT_DEATH_DATE"

## Column name that contains the date of event
if(!exists("death.col.name"))
    death.col.name <- "Date of death"

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


## The following code was necessary for the first time on 20200324. Hopefully it can be commented out and ignored in future iterations
dth.dat <- read_csv(paste0("../../../Data/Deaths/", date.data, " COVID19 Deaths.csv"))
dth.dat$dod <- lubridate::as_date(dth.dat$dod, format = "%d/%m/%Y", tz = "GMT")
dth.dat$dod_NHSE <- lubridate::as_date(dth.dat$dod_NHSE, format = "%m/%d/%Y", tz = "GMT")
x <- dth.dat$dod
x[is.na(x)] <- dth.dat$dod_NHSE[is.na(x)]
dth.dat <- dth.dat %>%
    mutate(Date = x)
## ## 

latest.date <- lubridate::as_date(date.data) - reporting.delay
earliest.date <- lubridate::as_date("2020-02-17")
all.dates <- as.character(seq(earliest.date, latest.date, by = 1))

dth.dat <- dth.dat %>%
    ## mutate(Date = lubridate::as_date(apply(dth.dat, 1, function(x) as.Date(as.character(x[death.col.name]), format = "%d/%m/%Y")))) %>%
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
            file = paste0("../../data/deaths/deaths", date.data, "_", reg, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
}
