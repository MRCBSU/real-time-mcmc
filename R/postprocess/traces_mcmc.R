# This script creates the traceplots and mcmc.RData## Run from within output directory
load("tmp.RData")

out.dir <- getwd()

## Location of this script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

## Where are various directories?
file.loc <- dirname(dirname(dirname(thisFile())))
proj.dir <- file.loc
Rfile.loc <- file.path(file.loc, "R/output")

source(file.path(Rfile.loc, "tracePlots.R"))
setwd(startwd)
