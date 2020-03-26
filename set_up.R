require(rmarkdown)

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

source(file.path(proj.dir, "set_up_inputs.R"))
source(file.path(proj.dir, "set_up_pars.R"))

## Make the output directory if necessary
flg.createfile <- !file.exists(out.dir)
if(flg.createfile) system(paste("mkdir", out.dir))

## Get the population sizes
require(readr)
require(tidyr)
pop <- read_csv(build.data.filepath("", "popn2018_all.csv"))
pop.input <- NULL
for(reg in regions){
    pop.full <- pop[pop$Name %in% ons.regions[[reg]] & !is.na(pop$Name), -(1:3), drop = FALSE]
    pop.full <- apply(pop.full, 2, sum)
    if(age.grps == "All")
        pop.input <- c(pop.input, pop.full["All ages"])
    }
## Remove spaces from region name.
regions <- gsub(" ", "_", regions, fixed = TRUE)

plain_document <- output_format(
    knitr = knitr_options(),
    pandoc = pandoc_options(to = "plain", ext = ".txt"),
)

pars.template.loc <- file.path(proj.dir, "inputs", "mod_pars.Rmd")
inputs.template.loc <- file.path(proj.dir, "inputs", "mod_inputs.Rmd")
render(pars.template.loc, output_dir = out.dir, output_format = plain_document)
render(inputs.template.loc, output_dir = out.dir, output_format = plain_document)

setwd(out.dir)
system("./rtm_gnu > runtime.txt", intern = TRUE)
