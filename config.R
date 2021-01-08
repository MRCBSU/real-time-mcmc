#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

# Either ONS or NHS
region.type <- "ONS"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(0)) %>% format("%Y%m%d"))
if (length(args) < 3) args <- c(args, "All", "England")

if (!exists("date.data")) date.data <- args[1]
stopifnot(!is.na(ymd(date.data)))
if (args[2] == "All")  {
	if (region.type == "NHS") 
		regions <- c("East_of_England", "London", "Midlands",
                     "North_East_and_Yorkshire", "North_West",
                     "South_East", "South_West"## ,
                     ## "Scotland", "Northern_Ireland", "Wales"
                     )
	else if (region.type == "ONS") 
		regions <- c("East_of_England", "East_Midlands",
					 "London", 
                     "North_East", "North_West",
                     "South_East", "South_West",
					 "West_Midlands", "Yorkshire_and_The_Humber"
                     )
	else stop("Unkown region type")
	nr <- length(regions)
} else {
	nr <- as.integer(args[2])
	regions <- args[3:(nr+2)]
	stopifnot(length(regions) == nr)
}

use.previous.run.for.start <- TRUE
previous.run.to.use <- "/home/jbb50/rds/hpc-work/real-time-mcmc/model_runs/20201127/NoPrev_varySero_cm4_ifr_60cutoff_matrices_20201127_deaths/"
iteration.number.to.start.from <- 6400

serology.delay <- 25 ## Assumed number of days between infection and developing the antibody response
sero.end.date <- ymd(20200605)

google.data.date <- format(ymd("20210108"), format = "%Y%m%d")
matrix.suffix <- "_timeuse_household_new_base"

## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
start.date <- lubridate::as_date("20200217")
nforecast.weeks <- 3
ndays <- as.integer(ymd(date.data) - start.date + (7 * nforecast.weeks) + 1)

cm.breaks <- seq(from = 36, to = ndays, by = 7) ## Day numbers where breaks happen

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
# adjusted_median: reporting-delay adjusted deaths produced by Pantelis, using medians
# adjusted_mean: reporting-delay adjusted deaths produced by Pantelis, using means
data.desc <- "deaths"

## The 'gp' stream in the code is linked to the pillar testing data
gp.flag <- 0	# 0 = off, 1 = on
## The 'hosp' stream in the code is linked to death data
hosp.flag <- 1					# 0 = off, 1 = on
## Do we want to include prevalence estimates from community surveys in the model?
prev.flag <- 0
prev.prior <- "Cevik" # "relax" or "long_positive" or "tight
## Shall we fix the serological testing specificity and sensitivty?
fix.sero.test.spec.sens <- FALSE #prev.flag == 1
exclude.eldest.prev <- FALSE

## Give the run a name to identify the configuratio
if (prev.flag) scenario.name <- paste0("Prev", prev.prior)
if (!prev.flag) scenario.name <- "NoPrev"
if (fix.sero.test.spec.sens) scenario.name <- paste0(scenario.name, "_fixedSero")
if (exclude.eldest.prev) scenario.name <- paste0(scenario.name, "_exclude_elderly_prev")

## Give the run a name to identify the configuratio
contact.model <- 4
if (contact.model != 4) scenario.name <- paste0(scenario.name, "_cm", contact.model) ## _latestart" ## _morefreq"
## Does each age group have a single IFR or one that varies over time?
single.ifr <- FALSE
if(single.ifr) scenario.name <- paste0(scenario.name, "_constant_ifr")

flg.confirmed <- (data.desc != "all")
flg.cutoff <- TRUE
if(flg.cutoff) {
	str.cutoff <- "60"
	scenario.name <- paste0(scenario.name, "_", str.cutoff, "cutoff")
}
if (data.desc == "all") {
	reporting.delay <- 18
} else if (data.desc == "reports") {
	reporting.delay <- 0
} else if (data.desc == "deaths") {
    reporting.delay <- 6
} else if (grepl("adjusted", data.desc)) {
    date.adj.data <- ymd(date.data)-2
    reporting.delay <- 1
} else {
	stop("Unknown data description")
}

data.dirs <- file.path(proj.dir,
                       paste0("data/RTM_format/", region.type, "/", c("deaths","serology","cases","prevalence"))                       
                       )
names(data.dirs) <- c("deaths", "sero", "cases", "prev")
      
flg.confirmed = TRUE

English.regions <- c("East_of_England", "London", "Midlands",
								  "North_East_and_Yorkshire", "North_West",
								  "South_East", "South_West", "England")
running.England <- any(regions %in% English.regions)


if(gp.flag){
    case.positivity <- TRUE ## include offsets for the number of tests taken
    ll.reporting.delay <- 4
    ## ll.start.date <- ymd("20200616") -- this is now defined in format_linelist.R
    ## Location where to find some incidence estimates immediately prior to the above-specified date
    outpp <- new.env()
    load(file.path(proj.dir, "model_runs", "20200619", "newContactModel6day_matrices_20200612_deaths", "output_matrices.RData"),
         envir = outpp)
    symptoms <- FALSE
    if(symptoms){
        scenario.name <- paste0(scenario.name, "_symptoms")
        asymptomatic.states <- "N"
    } else asymptomatic.states <- c("Y", "N", "U")
    pgp.prior.diffuse <- FALSE
} else case.positivity <- FALSE

## Get the date of the prevalence data
date.prev <- ymd("20201227")
## Convert that to an analysis day number
prev.end.day <- date.prev - start.date - 3
first.prev.day <- date.prev - 74 - start.date + 1
days.between.prev <- 14
## Default system for getting the days on which the likelihood will be calculated.
prev.lik.days <- rev(seq(from = as.integer(prev.end.day), to = as.integer(first.prev.day), by = -days.between.prev))
if(prev.flag) scenario.name <- paste0(scenario.name, "_prev", days.between.prev)

## ## Choose the name of the subdirectory in model_runs to use
out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
                     paste0(
							scenario.name,
							"_matrices_", google.data.date, matrix.suffix,
							"_", data.desc))	# Value actually used
if (!hosp.flag) out.dir <- paste0(out.dir, "_no_deaths")
if (gp.flag) out.dir <- paste0(out.dir, "_with_linelist")

threads.per.regions <- 2
