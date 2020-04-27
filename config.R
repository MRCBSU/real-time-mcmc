#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(1)) %>% format("%Y%m%d"))
if (length(args) < 3) args <- c(args, "1", "England")

if (!exists("date.data")) date.data <- args[1]
nr <- as.integer(args[2])
regions <- args[3:(nr+2)]
if (regions[1] == "All") regions <- c("East_of_England", "London", "Midlands",
									  "North_East_and_Yorkshire", "North_West",
									  "South_East", "South_West")
stopifnot(length(regions) == nr)

reporting.delay <- 5

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

## ## Choose the name of the subdirectory in model_runs to use
## subdir.name <- paste0(date.data, "regions_alone")
data.desc <- "reports" # Set to "reports" if running by reporting date
scenario.name <- "variable_relax_ifr_prior_delayAnne_confirmed"
out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
                     "age_stratified_1")	# Value actually used
data.dirs <- file.path(proj.dir,
                       "data/RTM_format/deaths")
## combined.dir <- file.path(proj.dir,
##                           "model_runs",
##                           date.data,
##                           paste0(nr, "regions_", data.desc, "_", scenario.name)) ## subdir.name, "_OVERALL_")	# Value actually used

## Do we want to consider only confirmed cases
flg.confirmed <- TRUE
=======
>>>>>>> origin/report
