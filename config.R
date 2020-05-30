#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################
library(lubridate)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c((today() - days(1)) %>% format("%Y%m%d"))
if (length(args) < 3) args <- c(args, "1", "Lombardy")

if (!exists("date.data")) date.data <- args[1]
nr <- as.integer(args[2])
regions <- args[3:(nr+2)]
stopifnot(length(regions) == nr)

reporting.delay <- 5
serology.delay <- 25 ## Assumed number of days between infection and developing the antibody response

scenario.name <- "no_m"

google.data.date <- format(ymd("20200522"), format = "%Y%m%d")
## Number of days to run the simulation for.
## Including lead-in time, analysis of data and short-term projection
start.date <- lubridate::as_date("20200101")
if (scenario.name == "no_m") {
	m.param.breakpoints <- NULL
} else {
	m.param.breakpoints <- ymd("20200307")
}
nforecast.weeks <- 3
ndays <- lubridate::as_date(date.data) - start.date + (7 * nforecast.weeks) + 1

## What age groupings are being used?
age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+") ## "All ages"
nA <- length(age.labs)

# The 'gp' stream in the code is linked to hospitalised cases
gp.flag <- 0					# 0 = off, 1 = on
## The 'hosp' stream in the code is linked to death data
hosp.flag <- 1					# 0 = off, 1 = on


out.dir <- file.path(proj.dir,
                     "model_runs",
                     date.data,
					 paste0("Lombardy_early_start", scenario.name)
					 )

