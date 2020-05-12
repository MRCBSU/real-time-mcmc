suppressMessages(library(lubridate))
suppressMessages(library(tidyverse))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

if(!exists("date.data"))
    date.data <- (today() - days(1)) %>% format("%Y%m%d")

## Where to find the data, if NULL use command line argument
if(!exists("sero.loc"))
    sero.loc <- paste0(date.sero, " week", week.sero, "_seroprev__modellers.csv")

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
    )

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
