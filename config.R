#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(1)) %>% format("%Y%m%d"))
if (length(args) < 3) args <- c(args, "All", "England")

if (!exists("date.data")) date.data <- args[1]
if (args[2] == "All")  {
	regions <- c("East_of_England", "London", "Midlands",
									  "North_East_and_Yorkshire", "North_West",
									  "South_East", "South_West")
	nr <- length(regions)
} else {
	nr <- as.integer(args[2])
	regions <- args[3:(nr+2)]
	stopifnot(length(regions) == nr)
}

reporting.delay <- 5
serology.delay <- 25 ## Assumed number of days between infection and developing the antibody response

google.data.date <- format(ymd("20200515"), format = "%Y%m%d")
## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
nforecast.weeks <- 3
ndays <- lubridate::as_date(date.data) - lubridate::as_date("20200216") + (7 * nforecast.weeks)

## What age groupings are being used?
age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+") ## "All ages"
nA <- length(age.labs)

if(!exists("regions")) regions <- "England"

region.code <- "Eng"

# Possible values:
# deaths: confirmed deaths only, by date of death
# reports: confirmed deaths only, by date of reporting
# all: all deaths, by date of death
data.desc <- "deaths" # Set to "reports" if running by reporting date
scenario.name <- "serology_varSens_slowantibodies"

flg.confirmed <- (data.desc != "all")
if (data.desc == "all") {
	reporting.delay <- 18
} else if (data.desc == "reports") {
	reporting.delay <- 0
} else if (data.desc == "deaths") {
	reporting.delay <- 5
} else {
	stop("Unknown data description")
}

## ## Choose the name of the subdirectory in model_runs to use
## subdir.name <- paste0(date.data, "regions_alone")
out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
                     paste0(scenario.name, "_matrices_", google.data.date)) ## Value actually used
data.dirs <- file.path(proj.dir,
                       c("data/RTM_format/deaths",
                         "data/RTM_format/serology")
                       )
names(data.dirs) <- c("deaths", "sero")
      
flg.confirmed = TRUE
