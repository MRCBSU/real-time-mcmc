require(dplyr)
require(tidyr)
require(tidyverse)
require(lubridate)

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
proj.dir <- file.loc
source(file.path(proj.dir, "R/output/results_api.R"))

create.spim.table <- function(data, name) {
  qprobs <- c(0.01, 0.05, 0.25, 0.75, 0.95, 0.99, 0.5)
  region.data <- data %>%
    get.aggregated.quantiles("region", qprobs) %>%
    bind_rows(data %>%
      get.aggregated.quantiles(NULL, qprobs) %>%
      mutate(region = "England")
    ) %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "q",
      values_from = value
    ) %>%
    rename(
      `1st centile` = q0.01,
      `5th centile` = q0.05,
      `25th centile` = q0.25,
      `75th centile` = q0.75,
      `95th centile` = q0.95,
      `99th centile` = q0.99,
      Value = q0.5,
      Geography = region
    ) %>%
    mutate(
      Group = "PHE/Cambridge",
      Model = "Regional/age",
      CreationDate = ymd(date.data),
      Scenario = "Forecast",
      ValueType = name
    )
}

tbl_inf <- create.spim.table(cum_infections, "infections_cum")
tbl_inf_inc <- create.spim.table(infections, "infections_inc")
tbl_deaths <- create.spim.table(noisy_deaths, "death_inc_line")

dir.string <- file.path(proj.dir, paste0("date_", date.data))
if(!file.exists(dir.string)) system(paste("mkdir", dir.string))

bind_rows(tbl_inf, tbl_deaths, tbl_inf_inc) %>%
    mutate(
		   `Creation Day` = day(CreationDate),
		   `Creation Month` = month(CreationDate),
		   `Creation Year` = year(CreationDate),
		   `Day of Value` = day(date),
		   `Month of Value` = month(date),
		   `Year of Value` = year(date),
		   Geography = str_replace_all(Geography, "_", " ")
  ) %>%
  select(-c(CreationDate, date)) %>%
  write.csv(file.path(dir.string, format(Sys.time(), "%X.csv")), row.names = FALSE)

