library(readr)
library(dplyr)
library(lubridate)

## Get case data
date.data <- "20200313"
ll.dat <- read_csv("~/stats_drive/COVID-19/LineList/Anon- 20200313 Linelist of UK cases (Cumulative) 0900h.csv")
## Get population data
cur.dir <- getwd()
setwd("/home/phe.gov.uk/paul.birrell/Documents/PHE/stats/Wuhan 2019 Coronavirus/Data/population")
source("get_popn.R")
setwd(cur.dir)
UK.pop <- sum(as.integer(gsub(",","",as.character(unlist(pop[pop$Name == "UNITED KINGDOM", -(1:4)])))))

age.grps <- c(0, 5, 15, 25, 45, seq(60, 80, by = 10))
age.labs <- c("<5yr", "5-14", "15-24", "25-44", "45-59", "60-69", "70-79", "80+")
## Currently latest date for age-independent modelling 12/3, latest date for age-dependent modelling 10/3
latest.date <- lubridate::as_date(date.data)
earliest.date <- lubridate::as_date("2020-03-02")

ll.dat <- ll.dat %>%
    filter(!is.na(as.integer(`Cumulative count`))) %>%
    mutate(Date = as.Date(`CMO Announce Date`, format = "%d/%m/%Y")) %>%
    filter(!is.na(Date)) %>%
    filter(Date <= latest.date) %>%
    filter(Date >= earliest.date) %>%
    mutate(fDate = factor(Date))
levels(ll.dat$fDate) <- c(levels(ll.dat$fDate), as.character(seq(earliest.date - 14, earliest.date, by = 1)))



## ll.dat should be the maximal dataset now.
## Use these to get the reduced datasets

## Whole UK data.
rtm.dat <- ll.dat %>%
    group_by(fDate, .drop = FALSE) %>%
    summarise(count = n()) ## %>%
    ## mutate(fDate = as.factor(Date))## , levels = seq(earliest.date, latest.date, by = 1)))

rtm.dat$fDate <- lubridate::as_date(rtm.dat$fDate)

rtm.dat <- arrange(rtm.dat, fDate)

rtm.denom <- data.frame(fDate = rtm.dat$fDate,
                        count = rep(UK.pop, nrow(rtm.dat))
                        )


## Write rtm.dat and rtm.denom to data files
write.table(rtm.dat,
            file = paste0("../data/Linelist/linelist", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
write.table(rtm.denom,
            file = paste0("../data/Linelist/ll_denom", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)

## ll.dat <- ll.dat %>%
##     mutate(Age_Grp = cut(Age,
##                          breaks = c(age.grps, Inf),
##                          right = FALSE))
                         
