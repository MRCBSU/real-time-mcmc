#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
args <- c((today() - days(1)) %>% format("%Y%m%d"),
          "1",
          "England")

if(length(args) > 0) date.data <- args[1]
if(length(args) > 1) nr <- as.integer(args[2])
if(length(args) >= (nr + 2)){
    regions <- NULL
    for(i in 1:nr) regions <- c(regions, args[i + 2])
} else { ## Not enough region names supplied
    stop("Not enough region names supplied")
}

if (!exists("date.data"))
    date.data <- (today() - days(1)) %>% format("%Y%m%d") ## Yesterday, by default
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
if (!exists("region.index")){
    region.index <- 1
} else 
    region.index <- which(names(all.regions) %in% regions)

region.code <- "Eng"

## ## Choose the name of the subdirectory in model_runs to use
## subdir.name <- paste0(date.data, "regions_alone")
data.desc <- "reports" # Set to "reports" if running by reporting date
scenario.name <- "variable_relax_ifr_prior_delayAnne_confirmed"
out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
                     paste0(nr, "regions_", region.code, data.desc, "_", scenario.name))	# Value actually used
data.dirs <- file.path(proj.dir,
                       "data/RTM_format/deaths")
combined.dir <- file.path(proj.dir,
                          "model_runs",
                          date.data,
                          paste0(nr, "regions_", data.desc, "_", scenario.name)) ## subdir.name, "_OVERALL_")	# Value actually used

## Do we want to consider only confirmed cases
flg.confirmed <- TRUE
