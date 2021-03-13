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
bak.out.dir <- out.dir <- getwd()
source(file.path(proj.dir, "R/output/results_api.R"))
out.dir <- bak.out.dir

create.spim.table <- function(data, name, by = NULL) {
  qprobs <- c(0.01, 0.05, 0.25, 0.75, 0.95, 0.99, 0.5)
  region.data <- data %>%
    get.aggregated.quantiles(c("region", by), qprobs) %>%
    bind_rows(data %>%
      get.aggregated.quantiles(by, qprobs) %>%
      mutate(region = "England")
    ) %>%
	mutate(
      Geography = region,
	  ValueType = name,
	) %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "quantile",
      values_from = value
    )
}


create.spim.table.noagg <- function(data, name) {
  qprobs <- seq(from = 0.05, to = 0.95, by = 0.05)
  region.data <- data %>%
    get.aggregated.quantiles("region", qprobs) %>%
    pivot_wider(
      names_from = "quantile",
      names_prefix = "quantile",
      values_from = value
    ) %>%
    mutate(
      Geography = region,
      ValueType = name,
    )
}

#fl.proj <- file.path(out.dir, "projections_midterm.RData")
spi.proj <- file.path(out.dir, "forSPI.RData")
#load(fl.proj)
load(spi.proj)
tbl_inf_inc <- create.spim.table(infections, "infections_inc")
dir.string <- file.path(proj.dir, paste0("phe-nowcasts/date_", date.data))
if(!file.exists(dir.string)) system(paste("mkdir", dir.string))
Rt.Eng <- get.infection.weighted.Rt(Rt, infections, which(parameter.to.outputs %in% iterations.for.Rt))
require(abind)
Rtx <- abind(Rt.Eng, Rt, along = 3)
names(dimnames(Rtx)) <- names(dimnames(Rt))
dimnames(Rtx)$region[1] <- "England"
tbl_R <- create.spim.table.noagg(Rtx, "R") ## For full history

bind_rows(tbl_inf_inc, tbl_R) %>%
  write.csv(file.path(out.dir, paste0("JBC_RTM_outputs", format(Sys.time(), "%Y%m%d_prev.csv"))), row.names = FALSE)

