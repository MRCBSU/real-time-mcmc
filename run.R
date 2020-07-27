library(rmarkdown)

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

## Which code is being considered
gp.flag <- 0
hosp.flag <- 1
sero.flag <- 0
viro.flag <- 0

## If these files don't already exits, make them
dir.data <- "data"
data.files <- paste0(data.dirs["deaths"], "/", data.desc, date.data, "_", regions, "_", nA, "ag", ifelse(flg.confirmed, "CONF", ""), ".txt")
names(data.files) <- regions
if(sero.flag){
  serosam.files <- paste0(data.dirs["sero"], "/", date.data, "_", regions, "_", nA, "ag_samples.txt")
  seropos.files <- paste0(data.dirs["sero"], "/", date.data, "_", regions, "_", nA, "ag_positives.txt")
} else {
  serosam.files <- seropos.files <- NULL
}
if(format.inputs){
  if(data.desc == "reports") {
	  source(file.path(proj.dir, "R/data/format_death_reports.R"))
  } else if (running.England) {
	  source(file.path(proj.dir, "R/data/format_deaths.R"))
  }
  if ("Scotland" %in% regions) {
	  source(file.path(proj.dir, "R/data/format_Scottish_deaths.R"))
  }
  if ("Northern_Ireland" %in% regions) {
	  source(file.path(proj.dir, "R/data/format_ni_deaths.R"))
  }
  if ("Wales" %in% regions) {
	  source(file.path(proj.dir, "R/data/format_wales_deaths.R"))
  }
  if(sero.flag){
	  source(file.path(proj.dir, "R/data/format_sero.R"))
  }
}

## Set up the model specification.
source(file.path(proj.dir, "set_up.R"))

## Compile the code
if(compile.code) {
    system("make rtm_optim")
}

## Run the code
startwd <- getwd()
setwd(out.dir)
save.image("tmp.RData")
if(run.code){
    system(file.path(proj.dir, "rtm_optim_8ag"), intern = TRUE)
	 system("chmod a-w coda* NNI* posterior* adaptive*")
}

## Post processing the results.
Rfile.loc <- file.path(file.loc, "R/output")

if(run.outputs){
    source(file.path(Rfile.loc, "tracePlots.R"))
	external = FALSE
	render(
		file.path(Rfile.loc, 'report-updated.Rmd'),
		html_document(pandoc_args = "--self-contained"),
		output_dir = out.dir,
		clean = FALSE, intermediates_dir = out.dir
	)
}

## Return back to initial directory
setwd(startwd)
