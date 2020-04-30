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

system(paste("mkdir -p", out.dir))

## do we need to do formatting?
format.inputs <- TRUE

## Will code need to be recompiled?
compile.code <- FALSE

## Do we want to actually run the code?
run.code <- TRUE

## Do we want to automatically run post-processing R code?
run.outputs <- TRUE

## Do the required data files exist?? If not, create them
data.files <- paste0(data.dirs, "/", data.desc, date.data, "_", regions, "_", nA, "ag", ifelse(flg.confirmed, "CONF", ""), ".txt")

## If these files don't already exits, make them
if(format.inputs && !all(file.exists(data.files))){
  dir.data <- "data"
  if(data.desc == "deaths")
	  source("R/data/format_deaths.R")
  if(data.desc == "reports")
	  source("R/data/format_death_reports.R")
}

## Set up the model specification.
source("set_up.R")

## Compile the code
if(compile.code) {
    system("make rtm_optim")
	system("chmod a-w coda* NNI* posterior* adaptive*")
}

## Run the code
startwd <- getwd()
setwd(out.dir)
if(run.code){
    system(file.path(proj.dir, "rtm_optim"), intern = TRUE)
} else save.image("tmp.RData")

## Post processing the results.
Rfile.loc <- file.path(file.loc, "R/output")

if(run.outputs){
    source(file.path(Rfile.loc, "tracePlots.R"))
	rmarkdown::render(
		file.path(Rfile.loc, 'report-updated.Rmd'),
		html_document(pandoc_args = "--self-contained"),
		output_dir = out.dir,
		clean = FALSE, intermediates_dir = out.dir
	)
}

## Return back to initial directory
setwd(startwd)
