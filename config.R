#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

# Either ONS or NHS
region.type <- "NHS"

args <- commandArgs(trailingOnly = TRUE)


if (length(args) == 0) args <- c((today() - days(0)) %>% format("%Y%m%d"))

#if (length(args) == 0) args <- c((today() - days(6)) %>% format("%Y%m%d"))#Paul's line

#if (length(args) == 0) args <- c((today() - days(3)) %>% format("%Y%m%d"))

#if (length(args) == 0) args <- c((today() - days(1)) %>% format("%Y%m%d"))

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

## Do we want to use NHSBT data (1) or RCGP data (0)
NHSBT.flag <- 1
## Do we want to use Roche N (0) or Roche S (1) data
RocheS.flag <- 0
## Assumed number of days between infection and developing the antibody response
serology.delay <- 25
## Last date for which serology is used
sero.end.date <- ymd(date.data) ## ymd(20200522) ## ymd(20210920)
## Last date for which first wave serology is used
sero.end.1stwv <- ymd(20200522)
## Format of dates used in the serology data
sero.date.fmt <- "%d%b%Y"
## Fix values at prior means?
fix.sero.test.spec.sens <- FALSE #prev.flag == 1


google.data.date <- format(ymd("2021-11-12"), format = "%Y%m%d")
#matrix.suffix <- "_timeuse_household_new_base"
matrix.suffix <- "_stable_household_new_base"



## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
start.date <- lubridate::as_date("20200217")
earliest.date <- ymd("2020-02-17")
nforecast.weeks <- 3
ndays <- as.integer(ymd(date.data) - start.date + (7 * nforecast.weeks) + 1)

cm.breaks <- seq(from = 36, to = ndays, by = 7) ## Day numbers where breaks happen
time.to.last.breakpoint <- 18 ## From the current date, when to insert the most recent beta breakpoint.
sdpar <- 100
break.window <- 2 ## How many WEEKS between breakpoints in the model for the transmission potential.

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

num.prev.days <- 557

## Shall we fix the serological testing specificity and sensitivty?
exclude.eldest.prev <- FALSE

use.INLA.prev <- TRUE

## Any inputs here for the vaccination data (or even if there is any)
vacc.flag <- 1
## Format used for dates in the vaccination file
vac.date.fmt <- "%Y-%m-%d"

## Give the run a name to identify the configuratio
if (prev.flag) scenario.name <- paste0("PrevINLA", num.prev.days)
if (!prev.flag) scenario.name <- "NoPrev"
if (fix.sero.test.spec.sens) scenario.name <- paste0(scenario.name, "_fixedSero")
scenario.name <- paste0(scenario.name, "Sero", ifelse(NHSBT.flag, "NHSBT", "RCGP"), "_", ifelse(sero.end.date == sero.end.1stwv, "1stwv", "All"))
if (exclude.eldest.prev) scenario.name <- paste0(scenario.name, "_exclude_elderly_prev")

## Give the run a name to identify the configuration
contact.model <- 6
contact.prior <- "ons"
## if (contact.model != 4)

    scenario.name <- paste0(scenario.name, "_cm", contact.model, contact.prior) ## _latestart" ## _morefreq"
## Does each age group have a single IFR or one that varies over time?
single.ifr <- FALSE
NHS28.alt.ifr.prior <- TRUE && region.type == "NHS"
if(single.ifr) scenario.name <- paste0(scenario.name, "_constant_ifr")
if(!single.ifr) ifr.mod <- "5bp"   ## 1bp = breakpoint over June, 2bp = breakpoint over June and October, lin.bp = breakpoint in June, linear increase from October onwards.
## tbreaks.interval <- 365.25 / 4

scenario.name <- paste0(scenario.name, "_IFR", ifr.mod, "")

flg.confirmed <- (data.desc != "all")
flg.cutoff <- TRUE
if(flg.cutoff) {
	str.cutoff <- "28"
	scenario.name <- paste0(scenario.name, "_", region.type, str.cutoff, "cutoff")
}
## Does each age group have a single IFR or one that varies over time?
single.ifr <- FALSE
NHS28.alt.ifr.prior <- (str.cutoff == "60") && (region.type == "NHS")
if(single.ifr) scenario.name <- paste0(scenario.name, "_constant_ifr")
if(!single.ifr) ifr.mod <- "5bp"   ## 1bp = breakpoint over June, 2bp = breakpoint over June and October, lin.bp = breakpoint in June, linear increase from October onwards.
## tbreaks.interval <- 365.25 / 4
scenario.name <- paste0(scenario.name, "_IFR", ifr.mod, "")
scenario.name <- paste0(scenario.name, "_", time.to.last.breakpoint, "wk", break.window)
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

## Is there a previous MCMC from which we can take some initial values?
use.previous.run.for.start <- TRUE
if(use.previous.run.for.start){
    if(region.type == "NHS"){
        if(str.cutoff == "60")

            previous.run.to.use <- file.path(proj.dir, "model_runs", "20211029", paste0(c("Prev542SeroNHSBT_All_NHS60cutoff_IFR5bp_18wk2_prev14-0PHE_matrices_20211029",
                                                                                          "Prev542SeroNHSBT_All_NHS60cutoff_IFR5bp_18wk2_prev14-0PHE_matrices_20211029"), matrix.suffix, "_", data.desc, c("_chain2", ""))
                                             )
        else previous.run.to.use <- file.path(proj.dir, "model_runs", "20211029", paste0(c("Prev542SeroNHSBT_All_NHS28cutoff_IFR5bp_18wk2_prev14-0PHE_matrices_20211029",
                                                                                           "Prev542SeroNHSBT_All_NHS28cutoff_IFR5bp_18wk2_prev14-0PHE_matrices_20211029"), matrix.suffix, "_", data.desc, c("_chain2", ""))
                                              )
    } else if(region.type == "ONS")
        previous.run.to.use <- file.path(proj.dir, "model_runs", "20211029", paste0(c("Prev542SeroNHSBT_All_ONS60cutoff_IFR5bp_18wk2_prev14-0PHE_matrices_20211029", # _stable_household_deaths_chain2",
                                                                                      "Prev542SeroNHSBT_All_ONS60cutoff_IFR5bp_18wk2_prev14-0PHE_matrices_20211029") # _stable_household_deaths")
                                                                                  , matrix.suffix, "_", data.desc, c("_chain2", ""))
                                         )
    
}
iteration.number.to.start.from <- 5000

## From where will the various datasets be sourced?
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
prev.cutoff.days <- 2
prev.days.to.lose <- 0
## Convert that to an analysis day number



date.prev <- lubridate::ymd("20211110")

prev.end.day <- date.prev - start.date - (prev.cutoff.days - 1) ## Last date in the dataset
last.prev.day <- prev.end.day - prev.days.to.lose ## Which is the last date that we will actually use in the likelihood?
first.prev.day <- prev.end.day - num.prev.days + 1
days.between.prev <- 14

## Default system for getting the days on which the likelihood will be calculated.
prev.lik.days <- rev(seq(from = as.integer(last.prev.day), to = as.integer(first.prev.day), by = -days.between.prev))
if(prev.flag) scenario.name <- paste0(scenario.name, "_prev", days.between.prev, "-", prev.days.to.lose)


## # Using 24 here means that each Friday an extra break will be added 3.5 weeks before the Friday in question
## lag.last.beta <- 24 - 7
## if (lag.last.beta != 24) scenario.name <- paste0(scenario.name, "_last_break_", lag.last.beta, "_days")



## if (matrix.suffix != "_timeuse_household_new_base") pasteo(scenario.name, "_", matrix.suffix)
efficacies <- "PHE" ## current values can be 'Nick', 'Jamie', or 'SPIM', 'PHE'.
vac.design <- "cumulative" ## currently values can be 'cumulative' or 'incident'.
scenario.name <- paste0(scenario.name, efficacies)

## ## temporary line - for adding ad hoc names to the scenario
## scenario.name <- paste0(scenario.name, "_manufacturer")

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

threads.per.regions <- 1

########### VACCINATION OPTIONS ###########
vacc.flag <- 1 ## Do we have any vaccination data




str.date.vacc <- "20211111" ## Optional: if not specified will take the most recent data file.

vacc.lag <- 21
vac.overwrite <- FALSE
if(vacc.flag){
    start.vac <- 301+vacc.lag ## Gives the day number of the first date for which we have vaccination data
    end.vac <- ndays ## Gives the most recent date for which we have vaccination data - or projected vaccination numbers
}
## How many vaccinations can we expect in the coming weeks
## - this is mostly set for the benefit of projections rather than model fitting.


future.n <- (c(0.16, 0.16, rep(0.16, 4), rep(0.16, 5)) * 10^6) * (55.98 / 66.65)




## Approximate data at which delta became dominant strain
delta.date <- ymd("20210510")
