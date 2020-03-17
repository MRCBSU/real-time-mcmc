library(assertr)
library(lubridate)
library(tidyverse)

read_JH_csv_file <- function(filename) {
    read_csv(filename) %>%                
        filter(`Province/State` == "United Kingdom") %>%
        verify(nrow(.) == 1) %>%                                      # Should only have the UK's one row left
        select(-c(`Province/State`, `Country/Region`, Lat, Long)) %>% # Remove all columns except the dates
        gather("day", "cumulative") %>%                               # Move from a single row to a row per day
        mutate(day = mdy(day)) %>%                                    # Convert to actual date types
        arrange(day) %>%                                              # Order by day, shouldn't be necessary but just in case!
        mutate(count = cumulative - lag(cumulative, default = 0)) %>% # Calculate daily increase
        assert(within_bounds(0, max(.$cumulative)), count) %>%        # Check daily changes plausible
        verify(count <= cumulative) %>%
        verify(sum(count) == slice(., nrow(.))$cumulative) %>%
        select(c(day, count))
}

read_JH_csv_file("time_series_19-covid-Deaths.csv") %>%
    write_tsv("RTM_deaths.txt", col_names = FALSE)

#read_JH_csv_file("time_series_19-covid-Confirmed.csv") %>%
#    write_tsv("RTM_confirmed_cases.txt", col_names = FALSE)
