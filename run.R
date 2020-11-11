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
format.inputs <- FALSE

## Will code need to be recompiled?
compile.code <- FALSE

## Do we want to actually run the code?
run.code <- FALSE

## Do we want to automatically run post-processing R code?
run.outputs <- FALSE

## Which code is being considered
if(!exists("gp.flag")) gp.flag <- 1
if(!exists("hosp.flag")) hosp.flag <- 1
if(!exists("sero.flag")) sero.flag <- 1
if(!exists("viro.flag")) viro.flag <- 0
if(!exists("prev.flag")) prev.flag <- 0

if (region.type == "NHS") {
	source(file.path(proj.dir, "R/data/get_NHS_pop.R"))
} else if (region.type == "ONS") {
	source(file.path(proj.dir, "R/data/get_ONS_pop.R"))
} else stop("Unknown region type for population")

## If these files don't already exits, make them
data.files <- paste0(data.dirs["deaths"], "/",
                     data.desc,
                     date.data, "_",
                     regions, "_",
                     nA, "ag",
                     ifelse(flg.confirmed, "CONF", ""),
                     reporting.delay, "delay")
if (exists("flg.cutoff")){
    if(flg.cutoff) 
	data.files <- paste0(data.files, "cutoff", str.cutoff)
}
data.files <- paste0(data.files, ".txt")
names(data.files) <- regions
if(sero.flag){
  serosam.files <- paste0(data.dirs["sero"], "/", date.data, "_", regions, "_", nA, "ag_samples.txt")
  seropos.files <- paste0(data.dirs["sero"], "/", date.data, "_", regions, "_", nA, "ag_positives.txt")
} else {
  serosam.files <- seropos.files <- NULL
}
if(gp.flag){
    cases.files <- paste0(data.dirs["cases"], "/", date.data, "_", regions, "_", nA, "_pillar_2_", ifelse(symptoms, "symptoms", "all"), ".txt")
    denoms.files <- paste0(data.dirs["cases"], "/", date.data, "_", regions, "_", nA, "_popdenom.txt")
} else {
    cases.files <- NULL
    denoms.files <- NULL
}
if(prev.flag){
    prev.mean.files <- paste0(data.dirs["prev"], "/", date.prev, "_", regions, "_ons_meanlogprev.txt")
    prev.sd.files <- paste0(data.dirs["prev"], "/", date.prev, "_", regions, "_ons_sdlogprev.txt")
} else {
    prev.mean.files <- NULL
    prev.sd.files <- NULL
}
if(format.inputs){
  if(data.desc == "reports") {
	  source(file.path(proj.dir, "R/data/format_death_reports.R"))
  } else if (grepl("adjusted", data.desc)) {
	  source(file.path(proj.dir, "R/data/format_adjusted_deaths.R"))
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
  if(gp.flag){
	  source(file.path(proj.dir, "R/data/format_linelist.R"))
  }
  if(prev.flag){
      source(file.path(proj.dir, "R", "data", "format_prev.R"))
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
if(exists("outpp"))
    rm(outpp)
save.image("tmp.RData")
if(run.code){
    system(file.path(proj.dir, "rtm_optim"), intern = TRUE)
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
