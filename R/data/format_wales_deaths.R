suppressMessages(library(lubridate))
suppressMessages(library(readxl))
suppressMessages(library(tidyverse))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

# Given a vector, return the first element where predicate is true
# Returns NULL if not true for any
first.where.true <- function(vec, predicate) {
	true.at <- which(predicate(vec))
	if (length(true.at) == 0) return(NULL)
	index.to.use <- min(true.at)
	return(vec[index.to.use])
}	

if(!exists("date.data"))
    date.data <- (today() - days(1)) %>% format("%Y%m%d")

reporting.delay <- 7

####################################################################
## BELOW THIS LINE SHOULD NOT NEED EDITING
####################################################################

## Location of this script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else if (.Platform$GUI == "RStudio" || Sys.getenv("RSTUDIO") == "1") {
                # We're in RStudio
                return(rstudioapi::getSourceEditorContext()$path)
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

## Where are various directories?
if(!exists("file.loc")){
    file.loc <- dirname(thisFile())
    proj.dir <- dirname(dirname(file.loc))
    dir.data <- file.path(proj.dir, "data")
	source(file.path(file.loc, "utils.R"))
	source(file.path(proj.dir, "config.R"))
}

latest.date <- ymd(date.data) - reporting.delay
earliest.date <- ymd("2020-02-17")

tbl_raw <- read_csv(
	"https://api.coronavirus.data.gov.uk/v1/data?filters=areaName=Wales&structure={%22areaName%22:%22areaName%22,%22date%22:%22date%22,%22newDeaths28DaysByDeathDate%22:%22newDeaths28DaysByDeathDate%22,%22newDeaths28DaysByPublishDate%22:%22newDeaths28DaysByPublishDate%22}&format=csv",
	col_types = cols(
		areaName = col_character(),
		date = col_date(format = ""),
		newDeaths28DaysByDeathDate = col_double(),
		newDeaths28DaysByPublishDate = col_double()
	)
)
file_date <- max(tbl_raw$date)
stopifnot(file_date == ymd(date.data))
stopifnot(all(tbl_raw$areaName == "Wales"))

rtm.dat <- tbl_raw %>%
  filter(date <= latest.date, date >= earliest.date) %>%
  select(Date = date, n = newDeaths28DaysByDeathDate) %>%
  replace_na(list(n = 0)) %>%
  arrange(Date)

stopifnot(nrow(rtm.dat) == latest.date - earliest.date + 1)

## Write rtm.dat to data file
print(paste(
	"Writing to",
	data.files["Wales"],
	"(", sum(rtm.dat$n), "total deaths,", nrow(rtm.dat), "rows.)"
))

if (nA == 1) {
  rtm.dat %>%
    select(Date, n) %>%
    write_tsv(
      data.files["Wales"],
      col_names = FALSE
    )
} else {
  rtm.dat %>%
    mutate(n0 = 0) %>%
    select(Date, n0, n) %>%
    write_tsv(
      data.files["Wales"],
      col_names = FALSE
    )
}
