library(assertr)
library(readr)
library(dplyr)
library(lubridate)

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

## Map our names for columns (LHS) to data column names (RHS)
col.names <- list(
	Lab_Report_Date = "Lab_Report_Date",
	Onsetdate = "Onsetdate"
)
## YYYYMMDD string, used in filenames and reporting lag
date.data <- "20200325"
## How long should the reporting lag be?
## Suggestion: overlay today's and yesterday's data
reporting.lag <- 2


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
file.loc <- dirname(thisFile())
proj.dir <- dirname(dirname(file.loc))
dir.data <- file.path(proj.dir, "data")
source(file.path(file.loc, "utils.R"))

## Which columns are we interested in?
ll.col.args <- list()
ll.col.args[[col.names[["Lab_Report_Date"]]]] <- col_character()
ll.col.args[[col.names[["Onsetdate"]]]] <- col_character()
ll.cols <- do.call(cols, ll.col.args)	# Calling with a list so use do.call

## Read the file and rename columns
ll.dat <- read_csv(
		build.data.filepath(subdir = "raw", date.data, " - Anonymised Line List.csv"),
		col_types = ll.cols
	) %>%
	rename(!!!col.names)

## Get population data
source(file.path(proj.dir, "R", "data", "get_popn.R"))
UK.pop <- sum(as.integer(gsub(",","",as.character(unlist(pop[pop$Name == "UNITED KINGDOM", -(1:4)])))))

## age.grps <- c(0, 5, 15, 25, 45, seq(60, 80, by = 10))
## age.labs <- c("<5yr", "5-14", "15-24", "25-44", "45-59", "60-69", "70-79", "80+")
## ## Currently latest date for age-independent modelling 12/3, latest date for age-dependent modelling 10/3
latest.date <- ymd(date.data)-reporting.lag
earliest.date <- ymd("2020-02-17")
all.dates <- as.character(seq(earliest.date, latest.date, by = 1))

ll.dat <- ll.dat %>%
    filter(!is.na(Lab_Report_Date))

ll.dat <- ll.dat %>%
    mutate(Date = fuzzy_date_parse(Lab_Report_Date)) %>%
	## Check plausibility (should flag parsing errors)
	verify(Date >= ymd("2020-01-01")) %>%
	verify(Date <= today()) %>%
	## Remore dates outside period of interest
    filter(Date <= latest.date) %>%
    filter(Date >= earliest.date) %>%
    mutate(fDate = factor(Date))
levels(ll.dat$fDate) <- c(levels(ll.dat$fDate), all.dates[!(all.dates %in% levels(ll.dat$fDate))])

rtm.dat <- ll.dat %>%
    group_by(fDate, .drop = FALSE) %>%
    summarise(count = n())

rtm.dat$fDate <- lubridate::as_date(rtm.dat$fDate)

rtm.dat <- arrange(rtm.dat, fDate)


rtm.denom <- data.frame(fDate = rtm.dat$fDate,
                        count = rep(UK.pop, nrow(rtm.dat))
                        )


## Write rtm.dat and rtm.denom to data files
write.table(rtm.dat,
            file = build.data.filepath(subdir = "RTM_format", "linelist", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
write.table(rtm.denom,
            file = build.data.filepath(subdir = "RTM_format", "ll_denom", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)

## This appears to be broken
if(FALSE & "Onsetdate" %in% names(ll.dat))
    {
        ## ll.dat <- ll.dat %>%
        ##     mutate(Age_Grp = cut(Age,
        ##                          breaks = c(age.grps, Inf),
        ##                          right = FALSE))
                         
        ons.dat <- ll.dat %>%
            filter(!is.na(Onsetdate))
        
        ons.dat <- ons.dat %>%
            mutate(ODate = fuzzy_date_parse(Onsetdate)) %>%
            mutate(Interval = Date - ODate) %>%
            mutate(Truncate = ymd(date.data)-3 - ODate)
        ons.dat$Interval[ons.dat$Interval < 0] <- 0
        ## truncated gamma distribution available in heavy package
        require(heavy)
        require(survival)
        
        ## but ptgamma in heavy package doesn't handle Inf correctly (see Chris Jackson's email 21/02/2020), so redefine it here
        ptgamma <- function(q, shape, scale=1, truncation=1, lower.tail=TRUE){
                                        # vectorise everything  
            n <- max(length(q), length(shape), length(scale), length(truncation))
            q <- rep(q, length=n)
            shape <- rep(shape, length=n)
            scale <- rep(scale, length=n)
            truncation <- rep(truncation, length=n)
            res <- numeric(n)
            res[q==Inf] <- 1
            sub <- q != Inf
            res[sub] <- heavy::ptgamma(q=q[sub], shape=shape[sub], scale=scale[sub], truncation=truncation[sub], lower.tail=lower.tail)
            res
        }
        
        ## create environment for custom distribution to pass to flexsurvreg
        custom.tgamma <- list(name="tgamma",
                              pars=c("shape","scale"),
                              location="scale",
                              transforms=c(log, log),
                              inv.transforms=c(exp, exp),
                              inits=function(t) { c(1.5, (2/3)*median(t)) })
        
        
        require(flexsurv)
        
        ons.dat <- ons.dat %>%
            mutate(report.min = Interval + 0.01 - ((1 - (Interval == 0)) * 0.51))
        
        Sdf <- Surv(time = ons.dat$report.min,
                    time2 = ons.dat$Interval + 0.5,
                    type = "interval2")
        (gamfit <- flexsurvreg(Sdf ~ 1, data = ons.dat, dist = custom.tgamma, aux = list(truncation = ons.dat$Truncate + 0.5)))

        gam.mean <- prod(exp(gamfit$coefficients))
        gam.sd <- sqrt(prod(exp(gamfit$coefficients[c(1,2,2)])))
    }
