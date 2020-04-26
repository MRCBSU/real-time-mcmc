suppressMessages(library(tidyverse))

## Location of this script
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

###############################################################

## Load files
if(!exists("proj.dir")){
  file.loc <- dirname(thisFile())
  proj.dir <- dirname(dirname(file.loc))
}
if (!exists("out.dir")) source(file.path(proj.dir, "config.R"))
if (!exists("infections")) {
  output.required <- file.path(out.dir, "output_matrices.RData")
  if (file.exists(output.required)) {
    load(output.required)
  } else {
    source(file.path(proj.dir, "R/output/tidy_output.R"))
  }
}

# TODO: read where we start/stop/thin from config files
parameter.iterations <- seq(from = 20000, to = 50000-1, by = 1)
outputs.iterations <- seq(from = 20000, to = 50000-1, by = 10)
parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)

################################################################

num.regions <- length(dimnames(infections)$region)
num.ages <- length(dimnames(infections)$age)
num.days <- length(dimnames(infections)$date)

parse.percentage <- function(percentage) {
  return(as.numeric(sub("%", "", percentage)) / 100)
}

get.aggregated.quantiles <- function(data, by, quantiles) {
  return(
    data %>%
      apply(c(by, "iteration", "date"), sum) %>%
      apply(c(by, "date"), quantile, probs = quantiles) %>%
      as.tbl_cube(met_name = "value") %>%
      as_tibble() %>%
      rename(quantile = Var1) %>%
      mutate(
        date = lubridate::as_date(date),
        quantile = parse.percentage(quantile)
      )
  )
}
