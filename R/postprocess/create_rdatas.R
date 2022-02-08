# Create all RData files needed for post-processing
load("tmp.RData")

out.dir <- getwd()

external <- FALSE

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

min.iteration <- as.integer(commandArgs(TRUE)[1])
if (is.na(min.iteration)) min.iteration <- burnin
print(glue::glue("Starting at iteration {min.iteration}, which is {min.iteration-burnin} iterations after the original burnin"))
stopifnot(min.iteration >= burnin)
source(file.path(Rfile.loc, "tidy_output.R"))
