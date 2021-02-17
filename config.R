#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

# Either ONS or NHS
region.type <- "ONS"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(1)) %>% format("%Y%m%d"))
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


serology.delay <- 25 ## Assumed number of days between infection and developing the antibody response
sero.end.date <- ymd(20200522)

google.data.date <- format(ymd("2021-02-12"), format = "%Y%m%d")
matrix.suffix <- "_timeuse_household_new_base"

## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
start.date <- lubridate::as_date("20200217")
nforecast.weeks <- 3
ndays <- as.integer(ymd(date.data) - start.date + (7 * nforecast.weeks) + 1)

cm.breaks <- seq(from = 36, to = ndays, by = 7) ## Day numbers where breaks happen
time.to.last.breakpoint <- 11 ## From the current date, when to insert the most recent beta breakpoint.

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
prev.flag <- 1
prev.prior <- "Cevik" # "relax" or "long_positive" or "tight
## Shall we fix the serological testing specificity and sensitivty?
fix.sero.test.spec.sens <- FALSE #prev.flag == 1
exclude.eldest.prev <- FALSE

## Any inputs here for the vaccination data (or even if there is any)
vacc.flag <- 1

## Give the run a name to identify the configuratio
if (prev.flag) scenario.name <- paste0("Prev", prev.prior)
if (!prev.flag) scenario.name <- "NoPrev"
if (fix.sero.test.spec.sens) scenario.name <- paste0(scenario.name, "_fixedSero")
if (exclude.eldest.prev) scenario.name <- paste0(scenario.name, "_exclude_elderly_prev")


## Give the run a name to identify the configuration
contact.model <- 4
if (contact.model != 4) scenario.name <- paste0(scenario.name, "_cm", contact.model) ## _latestart" ## _morefreq"
## Does each age group have a single IFR or one that varies over time?
single.ifr <- FALSE
if(single.ifr) scenario.name <- paste0(scenario.name, "_constant_ifr")
if(!single.ifr) ifr.mod <- "lin.bp"   ## 1bp = breakpoint over June, 2bp = breakpoint over June and October, lin.bp = breakpoint in June, linear increase from October onwards.
scenario.name <- paste0(scenario.name, "_IFR", ifr.mod)
flg.confirmed <- (data.desc != "all")
flg.cutoff <- TRUE
if(flg.cutoff) {
	str.cutoff <- "60"
	scenario.name <- paste0(scenario.name, "_", region.type, str.cutoff, "cutoff")
}
scenario.name <- paste0(scenario.name, "_", time.to.last.breakpoint)
if (data.desc == "all") {
	reporting.delay <- 18
} else if (data.desc == "reports") {
	reporting.delay <- 0
} else if (data.desc == "deaths") {
    reporting.delay <- 6
} else if (grepl("adjusted", data.desc)) {
    date.adj.data <- ymd(date.data) - 1  ## accounting for the fact that the raw and adjusted death files may have different dates on them.
    reporting.delay <- 1
} else {
	stop("Unknown data description")
}

data.dirs <- file.path(proj.dir,
                       paste0("data/RTM_format/", region.type, "/", c("deaths","serology","cases","prevalence","vaccination"))
                       )
names(data.dirs) <- c("deaths", "sero", "cases", "prev", "vacc")

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
num.prev.days <- 57
prev.cutoff.days <- 5
## Convert that to an analysis day number
date.prev <- lubridate::ymd("20210213") # Set this to last date in dataset
prev.end.day <- date.prev - start.date - prev.cutoff.days + 1
last.prev.day <- prev.end.day ## Which is the last date that we will actually use in the likelihood?
first.prev.day <- prev.end.day - num.prev.days + 1
days.between.prev <- 14

## Default system for getting the days on which the likelihood will be calculated.
prev.lik.days <- rev(seq(from = as.integer(last.prev.day), to = as.integer(first.prev.day), by = -days.between.prev))
if(prev.flag) scenario.name <- paste0(scenario.name, "_prev", days.between.prev)
if(prev.flag && prev.cutoff.days != 0) scenario.name <- paste0(scenario.name, "_skip", prev.cutoff.days)

# Using 24 here means that each Friday an extra break will be added 3.5 weeks before the Friday in question
lag.last.beta <- 24 - 7*2
if (lag.last.beta != 24) scenario.name <- paste0(scenario.name, "_last_break_", lag.last.beta, "_days")

if (matrix.suffix != "_timeuse_household_new_base") pasteo(scenario.name, "_", matrix.suffix)

## ## Choose the name of the subdirectory in model_runs to use
out.dir <- file.path(proj.dir,
                     "model_runs",
                     paste0(date.data, "-debug"),
                     paste0(
                         scenario.name,
                         "change_start_matrices_", google.data.date, matrix.suffix,
                         "_", data.desc))	# Value actually used
if (!hosp.flag) out.dir <- paste0(out.dir, "_no_deaths")
if (gp.flag) out.dir <- paste0(out.dir, "_with_linelist")

use.previous.run.for.start <- TRUE
previous.run.to.use <- "/home/jbb50/rds/hpc-work/real-time-mcmc/model_runs/20210212/PrevCevik_IFRlin.bp_ONS60cutoff_11_prev7_skip5_last_break_10_days_matrices_20210212_timeuse_household_new_base_deaths"
iteration.number.to.start.from <- 6400

threads.per.regions <- 2

########### VACCINATION OPTIONS ###########
vacc.flag <- 1 ## Do we have any vaccination data
str.date.vacc <- date.data ## Optional: if not specified will take the most recent data file.
vacc.rdata <- file.path(proj.dir, "data", "RTM_format", region.type, "vaccination", paste0(region.type, "vacc", str.date.vacc, ".RData"))
vacc.lag <- 21
vac.overwrite <- FALSE
if(vacc.flag){
    start.vac <- 301+vacc.lag ## Gives the day number of the first date for which we have vaccination data
    end.vac <- ndays ## Gives the most recent date for which we have vaccination data - or projected vaccination numbers
}

