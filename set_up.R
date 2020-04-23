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
## ## Load required functions for reading in data
source(file.path(proj.dir, "set_up_inputs.R"))
source(file.path(proj.dir, "set_up_pars.R"))

## Make the output directory if necessary
flg.createfile <- !file.exists(out.dir)
if(flg.createfile){
    tmp.dir <- out.dir
    while(!file.exists(dirname(tmp.dir[1])))
        tmp.dir <- c(dirname(tmp.dir[1]), tmp.dir)
    for(i in 1:length(tmp.dir))
        system(paste("mkdir -p", deparse(tmp.dir[i])))
}
## Change the hard-wiring of the number of age groups
header <- readLines("src/RTM_Header.h")
intHea <- grep("NUM_AGE_GROUPS", header)
header[intHea] <- paste0("#define NUM_AGE_GROUPS (", nA, ")")
write(header, file = "src/RTM_Header.h")

## Get the population sizes
require(readr)
require(tidyr)
require(dplyr)
load(build.data.filepath("population", "pop_nhs.RData"))
pop.input <- NULL
for(reg in regions){
    reg.nhs <- get.nhs.region(reg)
	if (reg == "Scotland" && age.labs[1] == "All") {
		pop.input <- c(pop.input, 5438100)
	} else {
		pop.full <- pop[pop$Name %in% nhs.regions[[get.nhs.region(reg)]] & !is.na(pop$Name), -(1:3), drop = FALSE]
		pop.full <- apply(pop.full, 2, sum)
		if(age.labs[1] == "All"){
                    pop.input <- c(pop.input, pop.full["All ages"])
                } else {
                    pdf <- data.frame(age = as.numeric(names(pop.full)[-1]), count = pop.full[-1])
                    pdf <- pdf %>%
                        mutate(age.grp = cut(pdf$age, age.agg, age.labs, right = FALSE, ordered_result = T)) %>%
                        group_by(age.grp) %>%
                        summarise(count = sum(count))
                    pop.input <- c(pop.input, pdf$count)
                }
        }
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
