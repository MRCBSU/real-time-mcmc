require(dplyr)
require(tidyr)
require(tidyverse)
require(lubridate)

suppressMessages(extract <- R.utils::extract)

mod.version.no <- 1.4
med.term.flag <- TRUE
nowcast.flag <- FALSE
mod.name <- ifelse(mod.version.no < 1.3, "Regional/age", ifelse(mod.version.no >= 1.4, "deaths/ons", "deaths and pillar2"))
## Get rid of any backslashes from, the model name
mod.fl.name <- gsub("/", "_", mod.name)

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
file.loc <- dirname(thisFile())
proj.dir <- dirname(dirname(file.loc))
if(!exists("out.dir")) out.dir <- getwd()
source(file.path(proj.dir, "R/output/results_api.R"))
out.dir <- getwd()
proj.dir <- dirname(dirname(dirname(out.dir)))
load(file.path(out.dir, "forSPI.RData"))
out.dir <- getwd()
proj.dir <- dirname(dirname(dirname(out.dir)))
projections.file <- "projections_midterm.RData"
scen.text <- "MTP"
save.text <- "MTP"
## projections.file <- "projections_R2.1.RData"
## scen.text <- "MTP R2.1"
## save.text <- "MTP_R_2.1"
mtp.filter.date <- lubridate::ymd("20210807") ## ymd(date.data)
dir.string <- file.path(proj.dir, paste0("spi-forecasts/date_", date.data))
if(!file.exists(dir.string)) system(paste("mkdir", dir.string))
nweeks.midterm <- 11

create.spim.table <- function(data, name, by = NULL) {
  qprobs <- seq(from = 0.05, to = 0.95, by = 0.05)
  region.data <- data %>%
    get.aggregated.quantiles(c("region", by), qprobs) %>%
    bind_rows(data %>%
      get.aggregated.quantiles(by, qprobs) %>%
      mutate(region = "England")
    ) %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "Quantile ",
      values_from = value
    ) %>%
    mutate(
      Group = "PHE/Cambridge",
      Model = mod.name,
      CreationDate = ymd(date.data),
      Version = mod.version.no,
      Geography = region,
      ValueType = name,
      Value = `Quantile 0.5`
    )
}

create.spim.table.noagg <- function(data, name) {
  qprobs <- seq(from = 0.05, to = 0.95, by = 0.05)
  region.data <- data %>%
    get.aggregated.quantiles("region", qprobs) %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "Quantile ",
      values_from = value
    ) %>%
    mutate(
      Group = "PHE/Cambridge",
      Model = mod.name,
      CreationDate = ymd(date.data),
      Version = mod.version.no,
      Geography = region,
      ValueType = name,
      Value = `Quantile 0.5`
    )
}

create.spim.summary <- function(data, name, q = NULL) {
    if(is.null(q)) q <- seq(from = 0.05, to = 0.95, by = 0.05)
    qdata <- data %>% apply("region", quantile, probs = q) %>%
        as.tbl_cube(met_name = "value") %>%
        as_tibble() %>%
        rename(quantile = Var1,
               Geography = region) %>%
        mutate(quantile = parse.percentage(quantile)) %>%
        pivot_wider(id_cols = Geography,
                    names_from = quantile,
                    names_prefix = "Quantile ",
                    values_from = value) %>%
        mutate(
            Group = "PHE/Cambridge",
            Model = mod.name,
            CreationDate = ymd(date.data),
            ValueType = name
        )
}

#### #### POPULATION MANIPULATION
## Aggregate youngest age groups in population table
population <- population %>%
    mutate(age2 = ifelse(age %in% c("<1yr", "1-4"), "0-4", age)) %>%
    group_by(region, age2) %>%
    summarise(popn = sum(value)) %>%
    rename(age = age2)

population <- population %>%
    bind_rows(population %>%
              group_by(age) %>%
              summarise(popn = sum(popn)) %>%
              mutate(region = "England"))

## Aggregate populations over age
population_reg <- population %>%
    group_by(region) %>%
    summarise(popn = sum(popn))
#### #### #### #### #### #### ####

if(nowcast.flag){
    ## tbl_inf <- create.spim.table(cum_infections, "infections_cum")
    ## tbl_deaths <- create.spim.table(noisy_deaths, "type28_death_inc_line")
    Rt.Eng <- get.infection.weighted.Rt(Rt, infections, which(parameter.to.outputs %in% iterations.for.Rt))
    require(abind)
    Rtx <- abind(Rt.Eng, Rt, along = 3)
    names(dimnames(Rtx)) <- names(dimnames(Rt))
    dimnames(Rtx)$region[1] <- "England"
        
    if(!exists("forecast.date")) forecast.date <- as_date("20201001")  ## today()-2
    forecast.idx <- which(dates.used >= forecast.date)  ## Want six weeks, so take seven to account for potentially mis-matching start date

    tbl_R <- create.spim.table.noagg(Rtx, "R") ## For full history
    ## tbl.R <- create.spim.summary(Rtx[, forecast.idx, ], "R") ## For current value
    NNI.tot <- apply(infections[,forecast.idx,, , drop = FALSE], c("date", "iteration", "region"), sum)
    xnames <- names(dimnames(NNI.tot))
    NNI.tot <- abind(NNI.tot, apply(NNI.tot, c("date", "iteration"), sum))
    names(dimnames(NNI.tot)) <- xnames
    dimnames(NNI.tot)$region[length(dimnames(NNI.tot)$region)] <- "England"
    tbl.nni <- create.spim.table.noagg(NNI.tot, "incidence")
    if(exists("prevalence")){
        prev.tot <- apply(prevalence[,forecast.idx,, , drop = FALSE], c("date", "iteration", "region"), sum)
        xnames <- names(dimnames(prev.tot))
        prev.tot <- abind(prev.tot, apply(prev.tot, c("date", "iteration"), sum))
        names(dimnames(prev.tot)) <- xnames
        dimnames(prev.tot)$region[length(dimnames(prev.tot)$region)] <- "England"
        tbl.prev <- create.spim.table.noagg(prev.tot, "prevalence") %>%
            inner_join(population_reg, by = c("Geography" = "region"))
        tbl.prev <- tbl.prev %>%
            mutate_at(vars(starts_with("Quantile")), `/`, (tbl.prev$popn / 100))  %>%
            ## mutate_at(Value, `/`, (tbl.prev$popn / 100)) %>%
            select(-popn)
    } else tbl.prev <- NULL
    
    ## Generation time
    Egt.arr <- do.call(cbind, Egt)
    Vargt.arr <- do.call(cbind, Vargt)
    
    names(dimnames(Egt.arr)) <- names(dimnames(Vargt.arr)) <- c("iteration", "region")
    tbl.egt <- create.spim.summary(Egt.arr, "mean_generation_time")
    tbl.vargt <- create.spim.summary(Vargt.arr, "kappa")
    
    today.int <- as.integer(ymd(date.data)) + 7
    start.int <- as.integer(ymd("20201201"))
    wanted.days <- as.character(start.int:today.int)
}

growth.prequantile <- function(data) {
    data %>%
        aperm(c("date", "region", "iteration", "age")) %>%
        `[`(wanted.days,,,,drop = FALSE) %>%
        apply(c("region", "iteration", "date"), sum) %>%
        as.tbl_cube(met_name = "value") %>%
        as_tibble() %>%
        mutate(
            date = lubridate::as_date(date)
        ) %>%
        arrange(date) %>%
        group_by(region, iteration) %>%
        summarise(growth = ((lead(value) - lag(value)) / (2 * value)), date = date) %>%
        filter(!is.na(growth)) %>%
        group_by(region, date)
    }
growth.prequantile.overall <- function(data) {
    data %>%
        aperm(c("date", "region", "iteration", "age")) %>%
        `[`(wanted.days,,,,drop = FALSE) %>%
        sum.all() %>%
        as.tbl_cube(met_name = "value") %>%
        as_tibble() %>%
        mutate(
            date = lubridate::as_date(date)
        ) %>%
        arrange(date) %>%
        group_by(iteration) %>%
        summarise(growth = ((lead(value) - lag(value)) / (2 * value)), date = date) %>%
        filter(!is.na(growth)) %>%
        group_by(date)
}
quantile.func <- function(data){
    data %>% summarise('Quantile 0.05' = quantile(growth, 0.05),
                       'Quantile 0.1' = quantile(growth, 0.1),
                       'Quantile 0.15' = quantile(growth, 0.15),
                       'Quantile 0.2' = quantile(growth, 0.2),
                       'Quantile 0.25' = quantile(growth, 0.25),
                       'Quantile 0.3' = quantile(growth, 0.3),
                       'Quantile 0.35' = quantile(growth, 0.35),
                       'Quantile 0.4' = quantile(growth, 0.4),
                       'Quantile 0.45' = quantile(growth, 0.45),
                       'Quantile 0.5' = quantile(growth, 0.5),
                       'Quantile 0.55' = quantile(growth, 0.55),
                       'Quantile 0.6' = quantile(growth, 0.6),
                       'Quantile 0.65' = quantile(growth, 0.65),
                       'Quantile 0.7' = quantile(growth, 0.7),
                       'Quantile 0.75' = quantile(growth, 0.75),
                       'Quantile 0.8' = quantile(growth, 0.8),
                       'Quantile 0.85' = quantile(growth, 0.85),
                       'Quantile 0.9' = quantile(growth, 0.9),
                       'Quantile 0.95' = quantile(growth, 0.95)
                       ) %>%
        mutate(Geography = str_replace_all(region, "_", " ")) %>%
        select(Geography, everything()) %>%
        ungroup() %>%
        select(-region) %>%
        mutate(
            Group = "PHE/Cambridge",
            Model = mod.name,
            CreationDate = ymd(date.data),
            ValueType = "growth_rate")
}

calc.growth.rate <- function(data) {
    data %>%
        growth.prequantile() %>%
        bind_rows(growth.prequantile.overall(data) %>% mutate(region = "England")) %>%
        quantile.func()
}

calc.infec.growth.rate <- function(R,gt,infecs) {
  ## To get a gt for England, just taking the mean of the other regions. In practice there is no regional heterogeneity. Would need more careful though if regional differences introduced.
  gt <- do.call(cbind, gt) %>% cbind(apply((.), 1, mean)) ## Growth rates usually the same across regions, so reduce dimensionality
  rt <- R %>%
     get.infection.weighted.Rt(infecs, which(parameter.to.outputs %in% iterations.for.Rt)) %>%
     abind(R, along = 3, use.dnns = TRUE) %>%
     as.tbl_cube(met_name = "value") %>%
     as_tibble() %>%
     mutate(region = ifelse(region == "", "England", region))
  dimnames(gt)[[1]] <- unique(rt$iteration)
  dimnames(gt)[[2]][dim(gt)[2]] <- "England"
  names(dimnames(gt)) <- c("iteration", "region")
  gt <- gt %>%
     as.tbl_cube(met_name = "genTime") %>%
     as_tibble()
  rt %>%
      left_join(gt) %>%
      mutate(growth = log(value) / genTime,
             date = lubridate::as_date(date)) %>%
      filter(date %in% wanted.days) %>%
      group_by(date, region) %>%
      quantile.func() ## %>%
    ##  summarise(Median = median(growth),
    ##           `95% CrI (lower)` = quantile(growth, 0.025),
    ##           `95% CrI (upper)` = quantile(growth, 0.975)
    ##           ) %>%
    ## arrange(match(region, c("England", regions))) %>%
    ## mutate(Region = str_replace_all(region, "_", " ")) %>%
    ## select(Region, everything()) %>%
    ## ungroup() %>%
    ## select(-region,-date)   
}

if(nowcast.flag){
    ## deaths.growth.rates <- calc.growth.rate(deaths)
    infecs.growth.rates <- calc.infec.growth.rate(Rt, Egt, infections)
    infecs.dbl.times <- infecs.growth.rates %>%
        mutate_at(
            vars(-c('Geography', 'date', 'Group', 'Model', 'CreationDate', 'ValueType')),
            ~ifelse(. < 0, NA, log(2)/(.))
        ) %>%
        rename('Quantile 0.05' = 'Quantile 0.95',
               'Quantile 0.1' = 'Quantile 0.9',
               'Quantile 0.15' = 'Quantile 0.85',
               'Quantile 0.2' = 'Quantile 0.8',
               'Quantile 0.25' = 'Quantile 0.75',
               'Quantile 0.3' = 'Quantile 0.7',
               'Quantile 0.35' = 'Quantile 0.65',
               'Quantile 0.4' = 'Quantile 0.6',
               'Quantile 0.45' = 'Quantile 0.55',
               'Quantile 0.5' = 'Quantile 0.5',
               'Quantile 0.55' = 'Quantile 0.45',
               'Quantile 0.6' = 'Quantile 0.4',
               'Quantile 0.65' = 'Quantile 0.35',
               'Quantile 0.7' = 'Quantile 0.3',
               'Quantile 0.75' = 'Quantile 0.25',
               'Quantile 0.8' = 'Quantile 0.2',
               'Quantile 0.85' = 'Quantile 0.15',
               'Quantile 0.9' = 'Quantile 0.1',
               'Quantile 0.95' = 'Quantile 0.05') %>%
        mutate(
            Group = "PHE/Cambridge",
            Model = mod.name,
            CreationDate = ymd(date.data),
            ValueType = "doubling_time")
    
    tbl.forecast <- bind_rows(## tbl_inf,
                              ## tbl_deaths,
                              tbl_R) %>%
    mutate(
        `Creation Day` = day(CreationDate),
        `Creation Month` = month(CreationDate),
        `Creation Year` = year(CreationDate),
        `Day of Value` = day(date),
        `Month of Value` = month(date),
        `Year of Value` = year(date),
        Geography = str_replace_all(Geography, "_", " "),
        `ModelType` = "Deaths",
        Scenario = "Nowcast",
        AgeBand = "All"
    ) %>%
        select(-c(CreationDate, date, region))
    
    tbl.nowcast <- bind_rows(tbl.nni, tbl.prev) %>%
        mutate(Version = mod.version.no,
            `Creation Day` = day(CreationDate),
            `Creation Month` = month(CreationDate),
            `Creation Year` = year(CreationDate),
            `Day of Value` = day(date),
            `Month of Value` = month(date),
            `Year of Value` = year(date),
            Geography = str_replace_all(Geography, "_", " "),
            `Value` = `Quantile 0.5`,
            `ModelType` = "Deaths",
            Scenario = "Nowcast",
            AgeBand = "All"
        ) %>%
        select(-CreationDate, -region, -date)

    tbl.current <- bind_rows(tbl.egt, tbl.vargt, infecs.growth.rates, infecs.dbl.times) %>%
        mutate(
            `Version` = mod.version.no,
            `Creation Day` = day(CreationDate),
            `Creation Month` = month(CreationDate),
            `Creation Year` = year(CreationDate),
            `Day of Value` = ifelse(ValueType %in% c("kappa", "mean_generation_time"),
                                    day(CreationDate),
                                    day(date)),
            `Month of Value` = ifelse(ValueType %in% c("kappa", "mean_generation_time"),
                                    month(CreationDate),
                                    month(date)),
            `Year of Value` = ifelse(ValueType %in% c("kappa", "mean_generation_time"),
                                    year(CreationDate),
                                    year(date)),
            Geography = str_replace_all(Geography, "_", " "),
            `Value` = `Quantile 0.5`,
            `ModelType` = "Deaths",
            Scenario = "Nowcast",
            AgeBand = "All"
        ) %>%
        select(-CreationDate)

    bind_rows(tbl.forecast, tbl.nowcast, tbl.current) %>%
        select(-date) %>%
        write.csv(file.path(dir.string, paste0("PHE_CAM_", mod.fl.name, format(Sys.time(), "%Y%m%d"), "_Nowcast.csv")), row.names = FALSE)
    ## bind_rows(tbl.forecast, tbl.nowcast, tbl.current %>% filter(ValueType == "growth_rate")) %>%
    ##     select(-date) %>%
    ##     mutate(Scenario = "Timeseries") %>%
    ##     write.csv(file.path(dir.string, paste0("PHE_CAM_", mod.fl.name, format(Sys.time(), "%Y%m%d"), "_TimeSeries.csv")), row.names = FALSE)
}

trim_forecast_days <- function(arrIn){
    nd <- 1:((nweeks.midterm + 1) * 7)
    nds <- as.integer(min(mtp.filter.day.no-1, end.hosp)) + nd
    arrOut <- extract(arrIn, indices = list(nds[nds <= length(dimnames(arrIn)$date)]), dims = "date", drop = FALSE)
    arrOut <- apply(arrOut, c("region", "date", "iteration"), function(x) c(x[1] + x[2], x[-(1:2)]))
    names(dimnames(arrOut))[1] <- "age"
    arrOut
    }

### MEDIUM-TERM FORECASTING ###
if(med.term.flag){
    mtp.filter.day.no <- mtp.filter.date - start.date + 1
    fl.proj <- file.path(out.dir, projections.file)
    if(!file.exists(fl.proj))
        stop("Missing projections file")
    load(fl.proj)
    if(exists("prevalence"))
        if(all(dim(prevalence) == dim(infections)))
            prev.flag <- TRUE
    deaths <- trim_forecast_days(deaths)
    infections <- trim_forecast_days(infections)
    tbl_proj <- create.spim.table(infections, "infections_inc", by = "age") %>% mutate(age = recode(age, '<1yr' = "0-4"))
    tbl_dproj <- create.spim.table(deaths, "type28_death_inc_line", by = "age") %>% mutate(age = recode(age, '<1yr' = "0-4"))
    if(prev.flag){
        prevalence <- trim_forecast_days(prevalence)
        tbl_aproj <- create.spim.table(prevalence, "prevalence_mtp", by = "age") %>% mutate(age = recode(age, '<1yr' = "0-4")) %>%
            left_join(population)
        tbl_aproj <- tbl_aproj %>%
            mutate_at(vars(starts_with("Quantile")), `/`, (tbl_aproj$popn / 100)) %>%
            mutate_at(vars("Value"), `/`, (tbl_aproj$popn / 100)) %>%
            select(-popn)
    } else tbl_aproj <- NULL
    
    tbl_midterm_forecast <- bind_rows(tbl_proj, tbl_dproj, tbl_aproj) %>%
        filter(date >= mtp.filter.date) %>%
        mutate(
            `ModelType` = "Deaths",
            `Creation Day` = day(CreationDate),
            `Creation Month` = month(CreationDate),
            `Creation Year` = year(CreationDate),
            `Day of Value` = day(date),
            `Month of Value` = month(date),
            `Year of Value` = year(date),
            AgeBand = ifelse(age == "<1yr", "0-4", age),
            ) %>%
        select(-c(CreationDate, date, region, age))

    ## write.csv(tbl_midterm_forecast,
    ##           file.path(dir.string, paste0("PHE_CAM_midterms", format(Sys.time(), "%Y%m%d.csv"))),
    ##           row.names = FALSE)
    
    tbl_all_proj <- create.spim.table(infections, "infections_inc")
    tbl_all_dproj <- create.spim.table(deaths, "type28_death_inc_line")
    if(prev.flag){
        tbl_all_aproj <- create.spim.table(prevalence, "prevalence_mtp") %>%
            left_join(population_reg)
        tbl_all_aproj <- tbl_all_aproj %>%
            mutate_at(vars(starts_with("Quantile")), `/`, (tbl_all_aproj$popn / 100)) %>%
            mutate_at(vars("Value"), `/`, (tbl_all_aproj$popn / 100)) %>%
            select(-popn)
    } else tbl_all_aproj <- NULL
    
    ## Joint table with population table
    ## Divide quantiles by population sizes
    
    tbl_midterm_all <- bind_rows(tbl_all_proj, tbl_all_dproj, tbl_all_aproj) %>%
        filter(date > ymd(date.data)) %>%
        mutate(
            `ModelType` = "Deaths",
            `Creation Day` = day(CreationDate),
            `Creation Month` = month(CreationDate),
            `Creation Year` = year(CreationDate),
            `Day of Value` = day(date),
            `Month of Value` = month(date),
            `Year of Value` = year(date),
            AgeBand = "All"
            ) %>%
        select(-c(CreationDate, date, region))

    tbl_midterm_forecast %>%
        bind_rows(tbl_midterm_all) %>%
        mutate(Scenario = scen.text,
               Geography = str_replace_all(Geography, "_", " ")) %>%
        write.csv(file.path(dir.string, paste0("PHE_CAM_", mod.fl.name, format(Sys.time(), "%Y%m%d"), "_", gsub(" ", "", save.text, fixed = TRUE), ".csv")),
                  row.names = FALSE)
        
}
### ####################### ###
