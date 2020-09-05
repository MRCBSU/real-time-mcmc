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
                     "South_East", "South_West"## ,
                     ## "Scotland", "Northern_Ireland", "Wales"
                     )
	nr <- length(regions)
} else {
	nr <- as.integer(args[2])
	regions <- args[3:(nr+2)]
	stopifnot(length(regions) == nr)
}

serology.delay <- 25 ## Assumed number of days between infection and developing the antibody response

google.data.date <- format(ymd("20200828"), format = "%Y%m%d")
## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
start.date <- lubridate::as_date("20200217")
nforecast.weeks <- 3
ndays <- lubridate::as_date(date.data) - start.date + (7 * nforecast.weeks) + 1

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
## scenario.name <- "base_varSens_pillar2_ifr"
scenario.name <- "base_varSens"
## scenario.name <- "base_varSens_pillar2"
## scenario.name <- "base_varSens_ifr"
contact.model <- 3

## The 'gp' stream in the code is linked to hospitalised cases
gp.flag <- 0					# 0 = off, 1 = on
## The 'hosp' stream in the code is linked to death data
hosp.flag <- 1					# 0 = off, 1 = on
## Does each age group have a single IFR or one that varies over time?
single.ifr <- TRUE

flg.confirmed <- (data.desc != "all")
if (data.desc == "all") {
	reporting.delay <- 18
} else if (data.desc == "reports") {
	reporting.delay <- 0
} else if (data.desc == "deaths") {
    flg.cutoff <- FALSE
    if(flg.cutoff) str.cutoff <- "60cod"
    reporting.delay <- 6
} else {
	stop("Unknown data description")
}
scenario.name <- paste0(scenario.name, reporting.delay, "day")


## ## Choose the name of the subdirectory in model_runs to use
## subdir.name <- paste0(date.data, "regions_alone")
out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
                     paste0(scenario.name, "_matrices_", google.data.date,
							"_", data.desc))	# Value actually used
if (!hosp.flag) out.dir <- paste0(out.dir, "_no_deaths")
if (gp.flag) out.dir <- paste0(out.dir, "_with_linelist")
data.dirs <- file.path(proj.dir,
                       c("data/RTM_format/deaths",
                         "data/RTM_format/serology",
                         "data/RTM_format/cases")
                       )
names(data.dirs) <- c("deaths", "sero", "cases")
      
flg.confirmed = TRUE

English.regions <- c("East_of_England", "London", "Midlands",
								  "North_East_and_Yorkshire", "North_West",
								  "South_East", "South_West", "England")
running.England <- any(regions %in% English.regions)


if(gp.flag){
    ll.reporting.delay <- 4
    ll.start.date <- ymd("20200616")
    ## Location where to find some incidence estimates immediately prior to the above-specified date
    outpp <- new.env()
    load(file.path(proj.dir, "model_runs", "20200619", "newContactModel6day_matrices_20200612_deaths", "output_matrices.RData"),
         envir = outpp)

}
