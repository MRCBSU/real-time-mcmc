library(tidyverse)
library(lubridate)
library(cubelyr)

snaps_dates <- ymd("20220729")

load("tmp.RData")
om <- new.env()

load("output_matrices.RData",envir = om)

## x <- state %>%
##     filter(state.names %in% c("S", "WS")) %>%
##     pivot_wider(id_cols = c("iteration", "region", "dose.numbers", "age.labs", "popn"), names_from = state.names, values_from = value)

df.ci <- om$cum_infections %>%
    as.tbl_cube(met_name = "cum_infections") %>%
    as_tibble() %>%
    inner_join(om$population) %>%
    rename(popn = value) %>%
    mutate(date = as_date(date)) %>%
    group_by(iteration, date) %>%
    summarise(cum_infections = sum(cum_infections)) %>%
    group_by(iteration) %>%
    summarise(infections = diff(cum_infections), date = date[-1])

df.ar <- om$AR %>%
    as.tbl_cube(met_name = "AR") %>%
    as_tibble() %>%
    inner_join(om$population) %>%
    rename(popn = value) %>%
    mutate(date = as_date(date)) %>%
    group_by(iteration, date) %>%
    summarise(cum_infected = sum(AR * popn)) %>%
    group_by(iteration) %>%
    summarise(infected = diff(cum_infected), date = date[-1])

df.overall <- inner_join(df.ci, df.ar) %>%
    filter(date %in% snaps_dates) %>%
    mutate(first_time_fraction = infected / infections) %>%
    group_by(date) %>%
    summarise(first_time_fraction = quantile(first_time_fraction, probs = c(0.025, 0.5, 0.975)))

## pct <- quantile(df.overall$first_time_pct, probs = c(0.025, 0.5, 0.975))

    
