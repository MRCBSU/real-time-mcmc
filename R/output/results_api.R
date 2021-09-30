library(coda)
library(cubelyr)
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
## Do we need to re-run all the calculations?
if (!exists("infections")) {
  if (!exists("out.dir")) {
	  warning("Importing config.R as no out.dir found")
	  source(file.path(proj.dir, "config.R"))
  }
  output.required <- file.path(out.dir, "output_matrices.RData")
  if (file.exists(output.required)) {
	mcmc_env <- tmp_env <- new.env()
    load(file.path(out.dir, "mcmc.RData"), mcmc_env)
    load(file.path(out.dir, "tmp.RData"), tmp_env)
    load(output.required)
    int_iter <- 0:(tmp_env$num.iterations - 1)
    parameter.iterations <- int_iter[(!((int_iter + 1 - tmp_env$burnin) %% tmp_env$thin.params)) & int_iter >= tmp_env$burnin]
    outputs.iterations <- int_iter[(!((int_iter + 1 - tmp_env$burnin) %% tmp_env$thin.outputs)) & int_iter >= tmp_env$burnin]
    parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)
    iterations.for.Rt <- parameter.to.outputs[seq(from = 1, to = length(parameter.to.outputs), length.out = 500)]
    stopifnot(length(parameter.to.outputs) == length(outputs.iterations)) # Needs to be subset
    date.data = mcmc_env$date.data
    dates.used= mcmc_env$dates.used
    start.date = mcmc_env$start.date
    end.hosp = mcmc_env$end.hosp	
        rm(int_iter)
	rm(tmp_env)
  } else {
    source(file.path(proj.dir, "R/output/tidy_output.R"))
  }
}

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
      add.quantiles(by, quantiles)
  )
}

get.infection.weighted.Rt <- function(R, infection, pars.to.out){
    inf_by_region <- apply(infection, c("iteration", "date", "region"), sum)[pars.to.out, , , drop = FALSE]
    Rt <- aperm(R, names(dimnames(inf_by_region)))
    stopifnot(all.equal(dim(Rt), dim(inf_by_region), check.names = FALSE))
    weighted_Rt_sum <- apply(inf_by_region * Rt, c("iteration", "date"), sum)
    inf_by_day <- apply(inf_by_region, c("iteration", "date"), sum)
    stopifnot(all.equal(dim(weighted_Rt_sum), dim(inf_by_day)))
    R.out <- weighted_Rt_sum / inf_by_day
    ## names(dimnames(R.out)) <- names(dimnames(R))
    R.out
}

add.quantiles <- function(data, by, quantiles) {
  return(    
    apply(data, c(by, "date"), quantile, probs = quantiles) %>%
      as.tbl_cube(met_name = "value") %>%
      as_tibble() %>%
      rename(quantile = Var1) %>%
      mutate(
        date = lubridate::as_date(date),
        quantile = parse.percentage(quantile)
      )
  )
}

matrix.to.tbl <- function(mat) {
  mat %>%
    as.tbl_cube(met_name = "value") %>%
    as_tibble() %>%
    mutate(
      date = lubridate::as_date(date),
      quantile = parse.percentage(quantile)
    )
}

sum.all <- function(arr) {
  return(apply(arr, c("iteration", "date"), sum))
}

sum.all.data <- function (df) {
  if (is.null(df)) return(NULL)
  df %>%
    group_by(date) %>%
    summarise(True = sum(value, na.rm = TRUE))
}
