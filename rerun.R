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

print("made it here")

if (!file.exists("mcmc.RData")) {
	source(file.path(Rfile.loc, "tracePlots.R"))
}

render(
	file.path(Rfile.loc, 'report-updated.Rmd'),
	html_document(pandoc_args = "--self-contained"),
	output_dir = out.dir
)

## Return back to initial directory
setwd(startwd)
