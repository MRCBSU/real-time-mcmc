#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(1)) %>% format("%Y%m%d"))
if (length(args) < 3) args <- c(args, 1, "Scotland")

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
English.regions <- c("East_of_England", "London", "Midlands",
								  "North_East_and_Yorkshire", "North_West",
								  "South_East", "South_West", "England")
running.England <- any(regions %in% English.regions)

serology.delay <- 25 ## Assumed number of days between infection and developing the antibody response

google.data.date <- format(ymd("20201120"), format = "%Y%m%d")
include.google <- TRUE
create.counterfactual <- FALSE
if (!create.counterfactual) {
  intervention.date <- NULL
} else {
  if (regions == "Wales") intervention.date <- ymd(20201030)
  if (regions == "Northern_Ireland") intervention.date <- ymd(20201013)
}

## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
start.date <- lubridate::as_date("20200217")
nforecast.weeks <- 9
ndays <- as.integer(lubridate::as_date(date.data) - start.date + (7 * nforecast.weeks) + 1)

## What age groupings are being used?
#age.labs <- c("All") ## "All ages"
age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+")
nA <- length(age.labs)
if (nA == 8) age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)

# Possible values:
# deaths: confirmed deaths only, by date of death
# reports: confirmed deaths only, by date of reporting
# all: all deaths, by date of death
data.desc <- "deaths"
# Possible values are NRS or PHS
scotland.data.desc <- "PHS"

scenario.name <- ifelse(nr == 1, regions, "")
scenario.name <- paste0(scenario.name, "_", nA, "ag", "_")
if (include.google) scenario.name <- paste0(scenario.name, "_with_google")
if (create.counterfactual) scenario.name <- paste0(scenario.name, "_with_intervention_end_date", intervention.date)
contact.model <- ifelse(nA == 1, 1, 3)

flg.confirmed <- (data.desc != "all")
if (data.desc == "all") {
	reporting.delay <- 18
} else if (data.desc == "reports") {
	reporting.delay <- 0
} else if (data.desc == "deaths") {
    flg.cutoff <- FALSE
    if(flg.cutoff) str.cutoff <- "60cod"
	if (running.England) {
      reporting.delay <- 14
	} else {
      if (regions == "Wales") reporting.delay <- 7
      if (regions == "Northern_Ireland") reporting.delay <- 2
      if (regions == "Scotland") reporting.delay <- 3
	}
} else {
	stop("Unknown data description")
}

# The 'gp' stream in the code is linked to hospitalised cases
gp.flag <- 0					# 0 = off, 1 = on
## The 'hosp' stream in the code is linked to death data
hosp.flag <- 1					# 0 = off, 1 = on


## ## Choose the name of the subdirectory in model_runs to use
## subdir.name <- paste0(date.data, "regions_alone")
out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
                     paste0(scenario.name, "_matrices_", google.data.date,
							"_", data.desc, "longer"))	# Value actually used
if (!hosp.flag) out.dir <- paste0(out.dir, "_no_deaths")
if (gp.flag) out.dir <- paste0(out.dir, "_with_hosp")
data.dirs <- file.path(proj.dir,
                       c("data/RTM_format/deaths",
                         "data/RTM_format/serology",
                         "data/RTM_format/cases")
                       )
names(data.dirs) <- c("deaths", "sero", "cases")
      
flg.confirmed = TRUE

