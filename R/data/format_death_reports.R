suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(lubridate))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

## YYYYMMDD string, used in filenames and reporting lag
# Default: yesterday's date
if(!exists("date.data"))
    date.data <- (today() - days(1)) %>% format("%Y%m%d")
# Or specify manually
# date.data <- "20200325"

## Where to find the data, if NULL use command line argument
if(!exists("deaths.loc"))
    deaths.loc <- paste0("Dataset Modelling " , date.data, ".csv") ## NULL

## Define an age-grouping
if(!exists("age.agg")){
    age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
    age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+")
    }
nA <- length(age.labs)

if(!exists("regions")){
    regions <- c("London", "Outside_London")
}

## Map our names for columns (LHS) to data column names (RHS)
col.names <- list(
    death_date = "dod",
    dbs_date = "DBSdeathreportdate",
    nhs_date = "NHSdeathreportdate",
    hpt_date = "HPTdeathreportdate",
    finalid = "finalid",
    onset_date = "onsetdate",
    nhs_region = "nhser_name",
    phe_region = "phec_name",
    utla_name = "utla_name"
)

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

## ## Where are various directories?
if(!exists("file.loc")){
    file.loc <- dirname(thisFile())
    proj.dir <- dirname(dirname(file.loc))
    dir.data <- file.path(proj.dir, "data")
    source(file.path(file.loc, "utils.R"))
}

if(!exists("data.files"))
    data.files <- build.data.filepath("RTM_format/deaths",
                                      "reports",
                                      date.data,
                                      "_",
                                      regions,
                                      "_",
                                      nA,
                                      "ages.txt")

## Which columns are we interested in?
death.col.args <- list()
death.col.args[[col.names[["death_date"]]]] <- col_character()
death.col.args[[col.names[["dbs_date"]]]] <- col_character()
death.col.args[[col.names[["nhs_date"]]]] <- col_character()
death.col.args[[col.names[["hpt_date"]]]] <- col_character()
death.col.args[[col.names[["onset_date"]]]] <- col_character()
death.col.args[[col.names[["finalid"]]]] <- col_double()
death.cols <- do.call(cols, death.col.args)	# Calling with a list so use do.call

plausible.death.date <- function(x) {
	within.range <- x$Report <= today() & x$Report >= ymd("2020-01-01")
        flgs <- is.na(x$Date)
        death.after.onset <- rep(TRUE, nrow(x))
        death.after.onset[flgs] <- (is.na(x$Onset) | (x$Onset <= x$Report))[flgs]
        death.after.onset[!flgs] <- (is.na(x$Onset) | (x$Onset <= x$Date))[!flgs]
        report.after.death <- is.na(x$Date) | (x$Date <= x$Report)
	return(within.range & death.after.onset & report.after.death)
}

## Some patients are known to have the month and day swapped
## Fix these here...
fix.dates <- function(df) {
	return(mutate(df,
		orig_date = Date,
		orig_onset = Onset,
                orig_report = Report,
		Date = heuristically.swap.day.and.month(Date),
		Onset = heuristically.swap.day.and.month(Onset),
                Report = heuristically.swap.day.and.month(Report)
	  ))
}

## Read the file and rename columns
if (is.null(deaths.loc)) {
	input.loc = commandArgs(trailingOnly = TRUE)[1]
} else {
	input.loc = build.data.filepath(subdir = "raw/deaths", deaths.loc)
}
print(paste("Reading from", input.loc))
dth.dat <- read_csv(input.loc,
                    col_types = death.cols) %>%
    rename(!!!col.names) %>%
    mutate(Date = fuzzy_date_parse(death_date),
           Onset = fuzzy_date_parse(onset_date),
           ReportDBS = fuzzy_date_parse(dbs_date),
           ReportNHS = fuzzy_date_parse(nhs_date),
           ReportHPT = fuzzy_date_parse(hpt_date)) %>%
    mutate(Report = pmin(ReportDBS, ReportNHS, ReportHPT, na.rm = TRUE)) %>%
    fix.dates %>%
    mutate(plausible_death_date = plausible.death.date(.) & !is.na(Report))


if (!all(dth.dat$plausible_death_date)) {
    implausible.dates <- dth.dat %>% filter(!plausible_death_date)
    print("WARNING: The following rows have implausible death dates and have been excluded: ")
    print(x.out <- implausible.dates %>% select(c(finalid, Onset, Date, Report)))## , swap_death, swap_onset)))
    dth.dat <- dth.dat %>% filter(plausible_death_date)
}

print("WARNING: the following rows have had onset and/or death date changed to become plausible: ")
dth.dat %>%
   filter(orig_date != Date | orig_onset != Onset | orig_report != Report) %>%
   select(c(finalid, Date, orig_date, Onset, orig_onset, Report, orig_report)) %>%
   mutate(
       orig_date = as_date(ifelse(orig_date == Date,
                                  NA, orig_date)),
       orig_onset = as_date(ifelse(orig_onset == Onset,
                                   NA, orig_onset)),
       orig_report = as_date(ifelse(orig_report == Report,
                                    NA, orig_report))
   ) %>%
   print(n=1000)

## The following code was necessary for the first time on 20200324. Hopefully it can be commented out and ignored in future iterations
#dth.dat$dod <- lubridate::as_date(dth.dat$dod, format = "%d/%m/%Y", tz = "GMT")
#dth.dat$dod_NHSE <- lubridate::as_date(dth.dat$dod_NHSE, format = "%m/%d/%Y", tz = "GMT")
#x <- dth.dat$dod
#x[is.na(x)] <- dth.dat$dod_NHSE[is.na(x)]
#dth.dat <- dth.dat %>%
    #mutate(Date = x)
## ## 

latest.date <- ymd(date.data)
earliest.date <- ymd("2020-02-17")
all.dates <- as.character(seq(earliest.date, latest.date, by = 1))

dth.dat <- dth.dat %>%
    filter(Report <= latest.date) %>%
    filter(Report >= earliest.date) %>%
    filter(age > 0) %>%  ## Large numbers of 0ys. Suspect age hasn't been recorded for the majority of these.
    ## mutate(fDate = factor(Report)) %>%
    mutate(Region = nhs.region(.)) %>%
    mutate(Region = map.to.region(Region)) %>%
    filter(Region %in% regions) %>%
    mutate(Age.Grp = cut(age, age.agg, age.labs, right = FALSE, ordered_result = T))
## levels(dth.dat$fDate) <- c(levels(dth.dat$fDate), all.dates[!(all.dates %in% levels(dth.dat$fDate))])

rtm.dat <- dth.dat %>%
    group_by(Report, Region, Age.Grp, .drop = FALSE) %>%
    tally %>%
    right_join(
        expand_grid(Report = as_date(earliest.date:latest.date), Region = regions, Age.Grp = age.labs),
        by = c("Report", "Region", "Age.Grp")
    ) %>%
    replace_na(list(n = 0)) %>%
    arrange(Report)

## rtm.dat$fDate <- lubridate::as_date(rtm.dat$fDate)

## rtm.dat <- rtm.dat %>%
##     arrange(fDate) %>%
##     filter(!is.na(Region))

## Write rtm.dat to data file
names(data.files) <- regions
for(reg in regions){
    region.dat <- pivot_wider(filter(rtm.dat, Region == reg),
                            id_cols = 1,
                            names_from = Age.Grp,
                            values_from = n)

    tmpFile <- data.files[reg]
    
    print(paste(
        "Writing to",
        tmpFile,
        "(", sum(region.dat[, -1]), "total deaths,", nrow(region.dat), "rows.)"))

    region.dat %>%
        write_tsv(
            tmpFile,
            col_names = FALSE
            )

}

## Save a quick plot of the data...
require(ggplot2)
rtm.dat %>%
    group_by(Report, Region) %>%
    summarise(count = sum(n)) -> rtm.dat.plot

gp <- ggplot(rtm.dat.plot, aes(x = Report, y = count, color = Region)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    ggtitle(paste("Daily number of death reports (on", lubridate::as_date(date.data), ")")) +
    xlab("Date of death report") +
    ylab("#Deaths") +
    theme(
        legend.position = "top",
        legend.justification = "left",
        )
ggsave(build.data.filepath("RTM_format/deaths", "reports_plot", date.data, ".pdf"),
       gp,
       width = 8.15,
       height = 6)
