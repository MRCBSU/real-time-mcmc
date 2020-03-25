library(readr)
library(dplyr)
library(lubridate)

## Get case data

## Inputs that should change on a daily basis
date.data <- "20200324"
dir.data <- "../../../Data/"

## Inputs that are dependent on the form of the data and the precise output required.
reporting.lag <- 2
## Get the format the dates are printed in in the input file.
date.format <- "%m/%d/%Y"

## Hopefully the following shouldn't change too frequently.
ll.dat <- read_csv(paste0(dir.data, "LineList/", date.data, " - Anonymised Line List.csv"))
## Get population data
cur.dir <- getwd()
setwd("../../../Data/population")
source("get_popn.R")
setwd(cur.dir)
UK.pop <- sum(as.integer(gsub(",","",as.character(unlist(pop[pop$Name == "UNITED KINGDOM", -(1:4)])))))

## age.grps <- c(0, 5, 15, 25, 45, seq(60, 80, by = 10))
## age.labs <- c("<5yr", "5-14", "15-24", "25-44", "45-59", "60-69", "70-79", "80+")
## ## Currently latest date for age-independent modelling 12/3, latest date for age-dependent modelling 10/3
latest.date <- lubridate::as_date(date.data)-reporting.lag
earliest.date <- lubridate::as_date("2020-02-17")
all.dates <- as.character(seq(earliest.date, latest.date, by = 1))

ll.dat <- ll.dat %>%
    filter(!is.na(Lab_Report_Date))

ll.dat <- ll.dat %>%
    mutate(Date = as.Date(Lab_Report_Date, format = date.format)) %>%
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
            file = paste0("../../data/Linelist/linelist", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)
write.table(rtm.denom,
            file = paste0("../../data/Linelist/ll_denom", date.data, ".txt"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)

if("Onsetdate" %in% names(ll.dat))
    {
        ## ll.dat <- ll.dat %>%
        ##     mutate(Age_Grp = cut(Age,
        ##                          breaks = c(age.grps, Inf),
        ##                          right = FALSE))
                         
        ons.dat <- ll.dat %>%
            filter(!is.na(Onsetdate))
        
        ons.dat <- ons.dat %>%
            mutate(ODate = lubridate::as_date(apply(ons.dat, 1, function(x) as.Date(as.character(x["Onsetdate"]), format = "%m/%d/%Y")))) %>%
            mutate(Interval = Date - ODate) %>%
            mutate(Truncate = (lubridate::as_date(date.data)-3) - ODate)
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
