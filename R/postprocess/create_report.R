## Run from within output directory
load("tmp.RData")

out.dir <- getwd()

require(rmarkdown)

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
file.loc <- dirname(thisFile())
proj.dir <- file.loc
Rfile.loc <- file.path(file.loc, "R/output")


# if (!file.exists("mcmc.RData")) {
	source(file.path(Rfile.loc, "tracePlots.R"))
        source(file.path(dirname(Rfile.loc), "postprocess", "create_rdatas.R"))
        print("created_rdatas")
# }


# suppressMessages(library(knitr))
# suppressMessages(library(DT))
# suppressMessages(library(lubridate))
# suppressMessages(library(plotly))
# suppressMessages(library(tidyverse))
# suppressMessages(extract <- R.utils::extract)
# cat(Rfile.loc)
# if(exists("Rfile.loc")) {
#   source(file.path(Rfile.loc, "results_api.R"))
#   source(file.path(Rfile.loc, "plot_funcs.R"))
# } else {
#   proj.dir <- "~/RTM"
#   source("results_api.R")
#   source("plot_funcs.R")
# }


render(
    file.path(Rfile.loc, 'report-updated.Rmd'),
    html_document(pandoc_args = "--self-contained"),
    output_dir = out.dir,
    intermediates_dir = "~/scratch/tmp"
)

# ## Return back to initial directory
# setwd(startwd)
