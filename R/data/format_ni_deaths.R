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


possible.files <- Sys.glob(glue::glue("/data/covid-19/data-raw/dstl/{ymd(date.data)}/{date.data}_All_SPIM_trust_*.xlsx"))
file.name <- possible.files[length(possible.files)]
rtm.dat <- read_excel(
  file.name,
  "Extracted Data",
  col_types = c("text", rep("numeric", 3), rep("text", 4), rep("numeric", 230)),
  na = c("", "*")
) %>%
  filter(Geography == "Northern Ireland") %>%
  transmute(
	Date = ymd(DateVal),
	n = SitRep_death_inc_line,
	Region = "Northern Ireland"
  ) %>%
  filter(Date <= latest.date, Date >= earliest.date) %>%
  replace_na(list(n = 0)) %>%
  arrange(Date)

stopifnot(nrow(rtm.dat) == latest.date - earliest.date + 1)

## Write rtm.dat to data file
print(paste(
	"Writing to",
	data.files["Northern_Ireland"],
	"(", sum(rtm.dat$n), "total deaths,", nrow(rtm.dat), "rows.)"
))

rtm.dat %>%
	select(Date, n) %>%
	write_tsv(
		data.files["Northern_Ireland"],
		col_names = FALSE
	)
