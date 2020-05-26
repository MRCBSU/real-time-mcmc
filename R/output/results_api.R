library(coda)
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
    load(file.path(out.dir, "mcmc.RData"))
    load(output.required)
    int_iter <- 0:(num.iterations - 1)
    parameter.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= burnin]
    outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]
    parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)
    iterations.for.Rt <- parameter.to.outputs[seq(from = 1, to = length(parameter.to.outputs), length.out = 500)]
    stopifnot(length(parameter.to.outputs) == length(outputs.iterations)) # Needs to be subset
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
    inf_by_region <- apply(infection, c("iteration", "date", "region"), sum)[pars.to.out, , ]
    Rt <- aperm(R, names(dimnames(inf_by_region)))
    stopifnot(all.equal(dim(Rt), dim(inf_by_region), check.names = FALSE))
    weighted_Rt_sum <- apply(inf_by_region * Rt, c("iteration", "date"), sum)
    inf_by_day <- apply(inf_by_region, c("iteration", "date"), sum)
    stopifnot(all.equal(dim(weighted_Rt_sum), dim(inf_by_day)))
    weighted_Rt_sum / inf_by_day
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
