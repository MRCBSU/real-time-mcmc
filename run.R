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
source(file.path(proj.dir, "config.R"))
source(file.path(proj.dir, "R/data/utils.R"))

## Will code need to be recompiled?
compile.code <- FALSE

## Do we want to actually run the code?
run.code <- FALSE

## Do we want to automatically run post-processing R code?
run.outputs <- FALSE

## Do the required data files exist?? If not, create them
data.files <- paste0(data.dirs, "/", data.desc, date.data, "_", regions, "ALL_", nA, "ag.txt")
## If these files don't already exits, make them
if(!all(file.exists(data.files))){
    dir.data <- "data"
    if(data.desc == "deaths")
        source("R/data/format_deaths.R")
    if(data.desc == "reports")
        source("R/data/format_death_reports.R")
}

## Set up the model specification.
source("set_up.R")

## Compile the code
if(compile.code)
    system("make rtm_optim")

## Run the code
startwd <- getwd()
setwd(out.dir)
if(run.code)
    system(file.path(proj.dir, "rtm_optim"), intern = TRUE)

## Post processing the results.
Rfile.loc <- file.path(file.loc, "R/output")

if(run.outputs){
    source(file.path(Rfile.loc, "tracePlots.R"))
    source(file.path(Rfile.loc, "projections.R"))
}

## Return back to initial directory
setwd(startwd)
