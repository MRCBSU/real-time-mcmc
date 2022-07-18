#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

# Either ONS or NHS
region.type <- "ONS"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(4)) %>% format("%Y%m%d"))
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


# Determine whether or not to run the model with the serology being dropped from a certain date onwards
# (Note: False -> doesn't drop the data)
Use_preprocessed_serology <- FALSE
preprocessed_sero_date <- ymd("20220627")

sero_cutoff_flag <- FALSE

if(sero_cutoff_flag) {
    #If dropping serology
    ## As date chosen, assumed to be chosen at a complete date
    serology.delay <- 25
    ## Last date for which serology is used
    sero.end.date <- ymd(20220301) ## ymd(20200522) ## ymd(20210920)

} else {
    #If not dropping serology
    ## Assumed number of days between infection and developing the antibody response
    serology.delay <- 25
    ## Last date for which serology is used
    sero.end.date <- ymd(date.data) ## ymd(20200522) ## ymd(20210920)
}
## Last date for which first wave serology is used
sero.end.1stwv <- ymd(20200522)
## Format of dates used in the serology data
sero.date.fmt <- "%d%b%Y"
## Fix values at prior means?
fix.sero.test.spec.sens <- FALSE #prev.flag == 1

# Flag to determine whether to cutoff the hospitalisation datastream early (T => use cutoff)
cutoff_hosps_early <- FALSE
# Variable to determine whether or not the admissions (T) or admissions + diagnoses (F) should be used
# Should nbe selected in combination with sus_seb_combination <- 3L in addition to having the preprocessed sus data
admissions_only.flag <- FALSE
## ## Value to note which combination of hospital data to use sus (0), sus + sebs (1), sebs only (2) or sus (preprocessed) + sebs (3)
sus_seb_combination <- 3L
## ##Value to note how many days to remove from the end of the dataset
adm_sus.strip_days <- 30L
adm_seb.strip_days <- 2L
seb_report_delay <- 1L  ## Used within this file, so can't be moved.
date.adm_seb <- ymd(20220709)
## date.adm_sus <- ymd(20210930)
date.adm.str <- lubridate::as_date(ifelse(sus_seb_combination > 0,
                                                  date.adm_seb - adm_seb.strip_days,
                                                  date.adm_sus - adm_sus.strip_days))
date_early_cutoff_hosps <- ymd(20220430)

## ## file.locs for admissions for geography linkers (with colname links)
adm.sus.geog_link.loc <- "utility_files/lad_to_region.csv"
adm.sus.geog_link <- "LAD19CD"
adm.sus.region_col <- "RGN19NM"
adm.seb.geog_link.loc <- "utility_files/trust lookup for paul.xlsx"
adm.seb.geog_link <- "Trust_code"
adm.seb.region_col <- "phec_nm"

## ## File names of pre-processed SUS data if it is to be used.
if(!admissions_only.flag) { ## Settings for preprocessed all hospitalisations
    adm_sus.end.date <- ymd(20201014) ## Date at which the sus and seb datasets are knitted together
    preprocessed_sus_names <- paste0("2022-01-02_", regions, "_6ag_counts.txt")
    sus_old_tab_sep <- T
    preprocessed_sus_csv_name <- "admissions_data_all_hosp.csv"
} else { ## Settings for admissions_only
    adm_sus.end.date <- ymd(20210505) ## Date at which the sus and seb datasets are knitted together
    preprocessed_sus_names <- paste0("2022-03-09_", regions, "_6ag_counts.txt")
    sus_old_tab_sep <- F
    preprocessed_sus_csv_name <- "admissions_data_admissions_only.csv"
}

names(preprocessed_sus_names) <- regions
print(preprocessed_sus_names)

## ## Admissions flags/dates
## adm.end.date <- date.data - adm_seb.strip_days ## Set this value if we want to truncate the data before its end.

## ## IF THE BELOW ARE NOT SPECIFIED, THE MOST RECENT FILE WILL BE CHOSEN
## date.adm_sus <- ymd()
## date.adm_seb <- ymd()

google.data.date <- format(ymd("20220708"), format = "%Y%m%d")
matrix.suffix <- "_timeuse_household"

## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
start.date <- lubridate::as_date("20200217")
earliest.date <- ymd("2020-02-17")
nforecast.weeks <- 3
ndays <- as.integer(ymd(date.data) - start.date + (7 * nforecast.weeks) + 1)

cm.breaks <- seq(from = 36, to = ndays, by = 7) ## Day numbers where breaks happen
time.to.last.breakpoint <- 11 ## From the current date, when to insert the most recent beta breakpoint.
sdpar <- 100
break.window <- 2 ## How many WEEKS between breakpoints in the model for the transmission potential.

## What age groupings are being used?
age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+") ## "All ages"
nA <- length(age.labs)

#! Added age groupings for the sitrep data
summarise_classes_sus <- list("[0,25)" = c("[0,1)", "[1,5)", "[5,15)", "[15,25)"))

summarise_classes_seb <- list("0_25" = c("0_5", "6_17", "18_24"),
                              "25_45"= c("25_34", "35_44"),
                              "45_65" = c("45_54", "55_64"),
                              "65_75" = c("65_74"),
                              "75" = c("75_84", "85")
)

age_adm.agg <- c(0, 25, 45, 65, 75, Inf) ## KEEP IN config.R
age_adm_seb.oldlabs <- c("0_25", "25_45", "45_65", "65_75", "75", "")
age_adm_sus.oldlabs <- c("[0,25)", "[25,45)", "[45,65)", "[65,75)", "[75,Inf]", "")
age_adm.labs <- c("0-25", "25-45","45-65", "65-75", "75+", "NA") ## "All ages"
nA_adm <- length(age_adm.labs)

#! Relabelling of regions for mismatches in the data (assumes consistency within dataset)
#! Additionally is run after get_region function therefore underscores are expected
relevel_sus_locs_nhs <- c("East_of_England" = "EAST_OF_ENGLAND", "London" = "LONDON", "Midlands" = "MIDLANDS",
                     "North_East_and_Yorkshire" = "NORTH_EAST_AND_YORKSHIRE", "North_West" = "NORTH_WEST",
                     "South_East" = "SOUTH_EAST", "South_West" = "SOUTH_WEST")

relevel_sus_locs_ons <- c()

relevel_seb_locs_ons <- c("Yorkshire_and_The_Humber" = "Yorkshire_and_Humber")

relevel_seb_locs_nhs <- c()

if(!exists("regions")) regions <- "England"

region.code <- "Eng"

## Possible values:
## deaths: confirmed deaths only, by date of death
## reports: confirmed deaths only, by date of reporting
## admissions: hospital admissions only, by date of admission
## all: all deaths, by date of death
## adjusted_median: reporting-delay adjusted deaths produced by Pantelis, using medians
## adjusted_mean: reporting-delay adjusted deaths produced by Pantelis, using means
data.desc <- "deaths"

## The 'gp' stream in the code is linked to the pillar testing data
gp.flag <- 0	# 0 = off, 1 = on
## Do we want the 'hosp' stream in the code linked to death data or to hospital admission data
deaths.flag <- hosp.flag <- 1			# 0 = admissions (by default - can be modified by explicitly setting adm.flag), 1 = deaths
## Do we want to include prevalence estimates from community surveys in the model?
prev.flag <- 1
prev.prior <- "Cevik" # "relax" or "long_positive" or "tight
num.prev.days <- 794
## Shall we fix the serological testing specificity and sensitivty?
exclude.eldest.prev <- FALSE

## Any inputs here for the vaccination data (or even if there is any)
vacc.flag <- 1
## Format used for dates in the vaccination file
vac.date.fmt <- "%d%b%Y"

## Deaths Flags
use_deaths_up_to_now_flag <- TRUE
custom_deaths_end_date <- lubridate::ymd("20220430")

Use_preprocessed_deaths <- FALSE

## Give the run a name to identify the configuratio
if (prev.flag) scenario.name <- paste0("Prev", num.prev.days)
if (!prev.flag) scenario.name <- "NoPrev"
if (fix.sero.test.spec.sens) scenario.name <- paste0(scenario.name, "_fixedSero")
scenario.name <- paste0(scenario.name, "Sero", ifelse(NHSBT.flag, "NHSBT", "RCGP"), "_", ifelse(sero.end.date == sero.end.1stwv, "1stwv", "All"))
if (exclude.eldest.prev) scenario.name <- paste0(scenario.name, "_exclude_elderly_prev")

## Give the run a name to identify the configuration
contact.model <- 6
contact.prior <- "ons"
## if (contact.model != 4)
    ## scenario.name <- paste0(scenario.name, "_cm", contact.model, contact.prior) ## _latestart" ## _morefreq"
flg.confirmed <- (data.desc != "all")
flg.cutoff <- TRUE
if(flg.cutoff) {
	str.cutoff <- ifelse(deaths.flag, ifelse(region.type == "ONS", "60", "28"), "")
	# str.cutoff <- "28"
	scenario.name <- paste0(scenario.name, "_", region.type, str.cutoff, "cutoff")
}
## Does each age group have a single IFR or one that varies over time?
single.ifr <- FALSE
NHS28.alt.ifr.prior <- (str.cutoff == "60") && (region.type == "NHS")
if(single.ifr) scenario.name <- paste0(scenario.name, "_constant_ifr")
if(!single.ifr) ifr.mod <- "8bp"   ## 1bp = breakpoint over June, 2bp = breakpoint over June and October, lin.bp = breakpoint in June, linear increase from October onwards.
## tbreaks.interval <- 365.25 / 4
scenario.name <- paste0(scenario.name, "_IFR", ifr.mod, "")
scenario.name <- paste0(scenario.name, ifelse(admissions_only.flag & data.desc == "admissions", "_admissions_only", ""), "_", time.to.last.breakpoint, "wk", break.window)
if (data.desc == "all") {
	reporting.delay <- 18
} else if (data.desc == "reports") {
	reporting.delay <- 0
} else if (data.desc == "deaths") {
    reporting.delay <- 6
} else if (grepl("adjusted", data.desc)) {
    date.adj.data <- ymd(date.data) - 1  ## accounting for the fact that the raw and adjusted death files may have different dates on them.
    reporting.delay <- 1
} else if (grepl("admissions", data.desc)) {
    reporting.delay  <- seb_report_delay
} else {
	stop("Unknown data description")
}

## Is there a previous MCMC from which we can take some initial values?
use.previous.run.for.start <- TRUE
if(use.previous.run.for.start){
    if(region.type == "NHS"){
        previous.run.to.use <- file.path(proj.dir, "model_runs", "20220701", paste0("Prev786SeroNHSBT_All_NHS", str.cutoff, "cutoff_IFR8bp_11wk2_prev14-0PHE_3dose_new_mprior_matrices_20220701", matrix.suffix, "_", ifelse(hosp.flag, "deaths", "admissions_no_deaths"), c("_chain2", ""))                                      )
    } else if(region.type == "ONS")
        previous.run.to.use <- file.path(proj.dir, "model_runs", ifelse(hosp.flag, "20220701", "20220702"), paste0("Prev786SeroNHSBT_All_ONS", str.cutoff, "cutoff_IFR8bp_", ifelse(admissions_only.flag & !hosp.flag, "admissions_only_", ""), "11wk2_prev14-0PHE_3dose_new_mprior_matrices_20220701", matrix.suffix, "_", ifelse(hosp.flag, "deaths", "admissions_no_deaths"), c("_chain2", ""))
                                         )
}

iteration.number.to.start.from <- 1

## From where will the various datasets be sourced?
#! Added admissions to data directories
data.dirs <- file.path(proj.dir,
                       paste0("data/RTM_format/", region.type, "/", c("deaths","serology","cases","prevalence","vaccination","admissions"))
                       )
names(data.dirs) <- c("deaths", "sero", "cases", "prev", "vacc", "adm")

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
prev.cutoff.days <- 4
prev.days.to.lose <- 0
## Convert that to an analysis day number
date.prev <- lubridate::ymd("20220705")
prev.end.day <- date.prev - start.date - (prev.cutoff.days - 1) ## Last date in the dataset
last.prev.day <- prev.end.day - prev.days.to.lose ## Which is the last date that we will actually use in the likelihood?
first.prev.day <- prev.end.day - num.prev.days + 1
days.between.prev <- 14

## Default system for getting the days on which the likelihood will be calculated.
prev.lik.days <- rev(seq(from = as.integer(last.prev.day), to = as.integer(first.prev.day), by = -days.between.prev))
if(prev.flag) scenario.name <- paste0(scenario.name, "_prev", days.between.prev, "-", prev.days.to.lose)

## if (matrix.suffix != "_timeuse_household_new_base") pasteo(scenario.name, "_", matrix.suffix)
efficacies <- "PHE" ## current values can be 'Nick', 'Jamie', or 'SPIM', 'PHE'.
vac.design <- "cumulative" ## currently values can be 'cumulative' or 'incident'.
scenario.name <- paste0(scenario.name, efficacies)

## ## temporary line - for adding ad hoc names to the scenario
## scenario.name <- paste0(scenario.name, "_manufacturer")

########### VACCINATION OPTIONS ###########
vacc.flag <- 1 ## Do we have any vaccination data
str.date.vacc <- "20220708" ## Optional: if not specified will take the most recent data file.
vacc.lag <- 21
vac.overwrite <- FALSE
if(vacc.flag){
    start.vac <- 301+vacc.lag ## Gives the day number of the first date for which we have vaccination data
    end.vac <- ndays ## Gives the most recent date for which we have vaccination data - or projected vaccination numbers
}
vac.n_doses <- 4L ## Number of doses in preprocessed data (Either 4, 3 or 2)
## How many vaccinations can we expect in the coming weeks
## - this is mostly set for the benefit of projections rather than model fitting.
future.n <- c(0.04, rep(0.04, 10)) * 10 ^ 6  * 55.98 / 66.65
future.booster.n <- c(1, 0.5, 0.3, 0.2, rep(0.2, 7)) * 10 ^ 6  * 55.98 / 66.65
future.fourth.n <- "c(1, 1, 1, 0.5, 0.5, rep(0.1, 6))* 10 ^ 6  * 55.98 / 66.65"
scenario.name <- paste0(scenario.name, "_", vac.n_doses, "dose")

## Approximate date at which delta became dominant strain (- one week)
delta.date <- ymd("20210503")
## Approximate date at which omicron became dominant strain (- one week)
omicron.date <- ymd("20211205")

scenario.name <- paste0(scenario.name, "_new_mprior")

## ## Choose the name of the subdirectory in model_runs to use
out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
                     paste0(
                         scenario.name,
                         ## Modified to rename the runs if cutting off the data early
                         ifelse(!use_deaths_up_to_now_flag & deaths.flag, paste0("_dropdeaths_", gsub("-", "",toString(custom_deaths_end_date))), ""),
                         ifelse(cutoff_hosps_early & !deaths.flag & !hosp.flag, paste0("_drophosp_", gsub("-", "",toString(date_early_cutoff_hosps))), ""),
                         ifelse(!sero_cutoff_flag, "", paste0("_dropsero_", gsub("-", "",toString(sero.end.date)))),
                         "_matrices2_", google.data.date, matrix.suffix,
                         "_", data.desc))	# Value actually used
if (!deaths.flag) out.dir <- paste0(out.dir, "_no_deaths")
if (gp.flag) out.dir <- paste0(out.dir, "_with_linelist")

threads.per.regions <- 1
