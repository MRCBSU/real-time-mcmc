#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

# Either ONS or NHS
region.type <- "NHS"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(1)) %>% format("%Y%m%d"))
if (length(args) < 3) args <- c(args, "All", "England")

if (!exists("date.data")) date.data <- args[1]
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
sero.end.date <- ymd(20200605)

google.data.date <- format(ymd("20201113"), format = "%Y%m%d")

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
# adjusted_median: reporting-delay adjusted deaths produced by Pantelis, using medians
# adjusted_mean: reporting-delay adjusted deaths produced by Pantelis, using means
data.desc <- "deaths"

## Give the run a name to identify the configuratio
scenario.name <- "Prev_relax_shortsero_fix_cm5"
contact.model <- 5

## The 'gp' stream in the code is linked to the pillar testing data
gp.flag <- 0	# 0 = off, 1 = on
## The 'hosp' stream in the code is linked to death data
hosp.flag <- 1					# 0 = off, 1 = on
## Do we want to include prevalence estimates from community surveys in the model?
prev.flag <- 1
## Does each age group have a single IFR or one that varies over time?
single.ifr <- FALSE
if(!single.ifr) scenario.name <- paste0(scenario.name, "_ifr")

if(!prev.flag) scenario.name <- paste0("No", scenario.name)
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
    date.adj.data <- ymd(date.data) - 1
    reporting.delay <- 1
} else {
	stop("Unknown data description")
}

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
        out.dir <- paste0(out.dir, "_symptoms")
        asymptomatic.states <- "N"
    } else asymptomatic.states <- c("Y", "N", "U")
    pgp.prior.diffuse <- FALSE
} else case.positivity <- FALSE

if(prev.flag){
    ## Get the date of the prevalence data
    date.prev <- ymd("20201109")
    date.old.prev <- ymd("20201028")
    date.start.prev <- ymd("20200802")
    ## Convert that to an analysis day number
    prev.end.day <- date.old.prev - start.date + 1 + 7
    prev.start.day <- date.start.prev - start.date
    ## Default system for getting the days on which the likelihood will be calculated.
    prev.lik.days <- rev(seq(as.integer(prev.end.day - 4), as.integer(prev.start.day), by = -28))
}

threads.per.regions <- 2
