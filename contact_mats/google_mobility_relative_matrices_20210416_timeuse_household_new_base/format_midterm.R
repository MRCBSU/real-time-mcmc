require(dplyr)
require(tidyr)
require(tidyverse)
require(lubridate)

mod.version.no <- 1.2
med.term.flag <- FALSE
mod.name <- ifelse(mod.version.no == 1.2, "Regional/age", "deaths and pillar2")

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
out.dir <- getwd()
proj.dir <- dirname(dirname(dirname(out.dir)))
source(file.path(proj.dir, "R/output/results_api.R"))
load(file.path(out.dir, "forSPI.RData"))
dir.string <- file.path(proj.dir, paste0("spi-forecasts/date_", date.data))
if(!file.exists(dir.string)) system(paste("mkdir -p", dir.string))

create.spim.table <- function(data, name) {
  qprobs <- seq(from = 0.05, to = 0.95, by = 0.05)
  region.data <- data %>%
    get.aggregated.quantiles("region", qprobs) %>%
	{ if (running.England) 
		bind_rows(., data %>%
					get.aggregated.quantiles(NULL, qprobs) %>%
					mutate(region = "England")
		)
	   else .
	} %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "Quantile ",
      values_from = value
    ) %>%
    mutate(
      Group = "PHE/Cambridge",
      Model = "Regional/age",
      CreationDate = ymd(date.data),
      Geography = region,
      ValueType = name,
      Value = `Quantile 0.5`
    )
}


fl.proj <- file.path(out.dir, "projections.RData")
if(!file.exists(fl.proj))
	stop("Missing projections file")
load(fl.proj)

tbl_all_proj <- create.spim.table(infections, "infections_inc")
tbl_all_dproj <- create.spim.table(deaths, "type28_death_inc_line")

tbl_midterm_all <- bind_rows(tbl_all_proj, tbl_all_dproj) %>%
	filter(date >= ymd(date.data)) %>%
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

output.file <- file.path(dir.string, paste0("PHE_RTM_midterms_", scenario.name, "_", format(Sys.time(), "%Y%m%d.csv")))
tbl_midterm_all %>%
	mutate(Geography = str_replace_all(Geography, "_", " ")) %>%
	write.csv(output.file, row.names = FALSE)
