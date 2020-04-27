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
  qprobs <- seq(from = 0.05, to = 0.95, by = 0.05)
  region.data <- data %>%
    get.aggregated.quantiles("region", qprobs) %>%
    bind_rows(data %>%
      get.aggregated.quantiles(NULL, qprobs) %>%
      mutate(region = "England")
    ) %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "Quantile ",
      values_from = value
    ) %>%
    mutate(
      Group = "PHE/Cambridge",
      Model = "Regional/age",
      CreationDate = ymd(date.data),
      Scenario = "Forecast",
      Geography = region,
      ValueType = name,
      Value = `Quantile 0.5`
    )
}

tbl_inf <- create.spim.table(cum_infections, "infections_cum")
tbl_deaths <- create.spim.table(noisy_deaths, "death_inc_line")

dir.string <- file.path(proj.dir, paste0("date_", date.data))
if(!file.exists(dir.string)) system(paste("mkdir", dir.string))

bind_rows(tbl_inf, tbl_deaths) %>%
    mutate(
		   `Creation Day` = day(CreationDate),
		   `Creation Month` = month(CreationDate),
		   `Creation Year` = year(CreationDate),
		   `Day of Value` = day(date),
		   `Month of Value` = month(date),
		   `Year of Value` = year(date),
		   Geography = str_replace_all(Geography, "_", " ")
  ) %>%
  select(-c(CreationDate, date, region)) %>%
  write.csv(file.path(dir.string, format(Sys.time(), "%X.csv")), row.names = FALSE)

