library(readr)
library(lubridate)
library(dplyr)

## Get deaths data
date.data <- "20200319"
dth.dat <- read_csv(paste0("../../../Data/Deaths/", date.data, " COVID19 Deaths.csv"))

latest.date <- lubridate::as_date(date.data) - 4
earliest.date <- lubridate::as_date("2020-02-17")
all.dates <- as.character(seq(earliest.date, latest.date, by = 1))

dth.dat <- dth.dat %>%
    mutate(Date = lubridate::as_date(apply(dth.dat, 1, function(x) as.Date(as.character(x["Date of death"]), format = "%d/%m/%Y")))) %>%
    filter(Date <= latest.date) %>%
    filter(Date >= earliest.date) %>%
    mutate(fDate = factor(Date)) %>%
    mutate(Region = as.factor(ifelse(Region == "London", "London", "Outside_London")))
levels(dth.dat$fDate) <- c(levels(dth.dat$fDate), all.dates[!(all.dates %in% levels(dth.dat$fDate))])

rtm.dat <- dth.dat %>%
    group_by(fDate, Region, .drop = FALSE) %>%
    summarise(count = n())

rtm.dat$fDate <- lubridate::as_date(rtm.dat$fDate)

rtm.dat <- arrange(rtm.dat, fDate)

## Write rtm.dat to data file
for(reg in levels(rtm.dat$Region)){
    write.table(filter(rtm.dat, Region == reg) %>%
               select(fDate, count),
            file = paste0("../../data/deaths/deaths", date.data, "_", reg, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
}
