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


## Will code need to be recompiled?
compile.code <- FALSE

## Do we want to actually run the code?
run.code <- TRUE

## Do we want to automatically run post-processing R code?
run.outputs <- TRUE

## Do the required data files exist?? If not, create them
data.files <- file.path(proj.dir, "data", "RTM_format", "Lombardy",
						paste0(date.data, ".txt"))

## Which code is being considered
gp.flag <- 0
hosp.flag <- 1
sero.flag <- 0
viro.flag <- 0


## Set up the model specification.
source(file.path(proj.dir, "set_up.R"))

## Compile the code
if(compile.code) {
    system("make rtm_optim")
}

## Run the code
startwd <- getwd()
setwd(out.dir)
if(run.code){
    system(file.path(proj.dir, "rtm_optim"), intern = TRUE)
	system("chmod a-w coda* NNI* posterior* adaptive*")
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
