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
proj.dir <- dirname(dirname(file.loc))
if(!exists("out.dir")) out.dir <- getwd()
out.dirx <- out.dir
source(file.path(proj.dir, "R/output/results_api.R"))



out.dir <- out.dirx; rm(out.dirx)


create.spim.table <- function(data, name, by = NULL) {
  qprobs <- c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975, 0.5)
  region.data <- data %>%
    get.aggregated.quantiles(c("region", by), qprobs) %>%
    bind_rows(data %>%
      get.aggregated.quantiles(by, qprobs) %>%
      mutate(region = "England")
    ) %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "q",
      values_from = value
    ) %>%
    rename(
      `1st centile` = q0.025,
      `5th centile` = q0.05,
      `25th centile` = q0.25,
      `75th centile` = q0.75,
      `95th centile` = q0.95,
      `99th centile` = q0.975,
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


fl.proj <- file.path(out.dir, "projections_midterm.RData")
load(fl.proj)
tbl_inf <- create.spim.table(cum_infections, "infections_cum")
tbl_inf_inc <- create.spim.table(infections, "infections_inc")
tbl_deaths <- create.spim.table(noisy_deaths, "death_inc_line")
tbl_deaths_age <- create.spim.table(deaths, "death_inc_line", by = "age")
tbl_inf_inc_age <- create.spim.table(infections, "infections_inc", by = "age")
## if(prev.flag){
    tbl_prev <- create.spim.table(prevalence, "prevalence")
    tbl_prev_age <- create.spim.table(prevalence, "prevalence", by = "age")
## } else {
    ## tbl_prev <- NULL
    ## tbl_prev_age <- NULL
## }
dir.string <- file.path(proj.dir, paste0("phe-nowcasts/date_", date.data))
if(!file.exists(dir.string)) system(paste("mkdir", dir.string))

bind_rows(tbl_inf, tbl_inf_inc, tbl_inf_inc_age, tbl_deaths, tbl_deaths_age, tbl_prev, tbl_prev_age) %>%
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
  write.csv(file.path(out.dir, paste0("PHE_RTM_outputs", format(Sys.time(), "%Y%m%d_prev.csv"))), row.names = FALSE)

