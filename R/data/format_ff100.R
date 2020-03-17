library(readr)
library(dplyr)
library(lubridate)

## Get case data
date.data <- "20200312"
ff.dat <- read_csv(paste0("~/Documents/PHE/stats/Wuhan 2019 Coronavirus/Epidemic_growth/UK/adjusted_cases_", date.data, ".csv"))
## Get population data
cur.dir <- getwd()
setwd("/home/phe.gov.uk/paul.birrell/Documents/PHE/stats/Wuhan 2019 Coronavirus/Data/population")
source("get_popn.R")
setwd(cur.dir)
UK.pop <- sum(as.integer(gsub(",","",as.character(unlist(pop[pop$Name == "UNITED KINGDOM", -(1:4)])))))

start.date <- lubridate::as_date("20200217")

## Transform Date to a factor with unused levels, to backfill.
ff.dat <- ff.dat %>%
    mutate(fDate = factor(Date))
levels(ff.dat$fDate) <- c(levels(ff.dat$fDate), as.character(seq(start.date, max(ff.dat$Date), by = 1)))

## Whole UK data for the RTM:
tempfun <- function(x){
    if(length(x) == 0) return(0)
    rem <- x[1] - floor(x[1])
    floor(x[1]) + ifelse(runif(1) < rem, 1, 0)
    }
rtm.dat <- ff.dat %>%
    group_by(fDate, .drop = FALSE) %>%
    summarise(count = tempfun(adjusted)) %>%
    arrange(lubridate::as_date(fDate))

## Trivial denominator - the full population are exposed
rtm.denom <- data.frame(fDate = rtm.dat$fDate,
                        count = rep(UK.pop, nrow(rtm.dat))
                        )

## Write rtm.dat and rtm.denom to data files
write.table(rtm.dat,
            file = paste0("../data/FF100/ff100_", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
write.table(rtm.denom,
            file = paste0("../data/FF100/ff_denom", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
