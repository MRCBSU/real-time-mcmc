library(tidyverse)
library(rmarkdown)
library(qs)

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

proj1 <- "projections_midterm.RData"
proj2 <- "projections_counter.RData"

mid.env <- new.env()
cou.env <- new.env()

load(proj1, envir = mid.env)

thin.outs <- function(data){
    data %>%
        R.utils::extract(indices = list(iters.keep), dims = "iteration") %>%
        apply(c("iteration", "age", "date"), sum)
}

iters.keep <- sample.int(length(dimnames(mid.env$deaths)[["iteration"]]), 1000)

deaths.vac <- thin.outs(mid.env$deaths)
infecs.vac <- thin.outs(mid.env$infections)
preval.vac <- thin.outs(mid.env$prevalence)
severe.vac <- thin.outs(mid.env$vacc.infections)

rm(mid.env)

load(proj2, envir = cou.env)

deaths.no.vac <- thin.outs(cou.env$deaths)
infecs.no.vac <- thin.outs(cou.env$infections)
preval.no.vac <- thin.outs(cou.env$prevalence)
severe.no.vac <- thin.outs(cou.env$vacc.infections)

Rfile.dir <- dirname(thisFile())
rm(cou.env)

## save(deaths.vac, deaths.no.vac, infecs.vac, infecs.no.vac, preval.vac, preval.no.vac, severe.vac, severe.no.vac, file = "deaths_comparison_prev.RData")
## estimate total deaths saved is 201 (164-251)!!

source(file.path(Rfile.dir, "vacc_lives_saved.R"))
## pth <- getwd()
## setwd(Rfile.dir)
## rmarkdown::render(
##     'lives_saved_report.Rmd',
##     output_dir = pth,
##     clean = T
## )
