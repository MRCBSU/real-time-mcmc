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
header.file <- file.path(proj.dir, "src/RTM_Header.h")
header <- readLines(header.file)
intHea <- grep("NUM_AGE_GROUPS", header)
header[intHea] <- paste0("#define NUM_AGE_GROUPS (", nA, ")")
write(header, file = header.file)

## Get the population sizes
pop.input <- c(75523, 335408, 962560, 947416, 2454341, 3012490, 1097720,
			   1175116)

regions <- gsub(" ", "_", regions, fixed = TRUE)

plain_document <- output_format(
    knitr = knitr_options(),
    pandoc = pandoc_options(to = "plain", ext = ".txt"),
)

pars.template.loc <- file.path(proj.dir, "inputs", "mod_pars.Rmd")
inputs.template.loc <- file.path(proj.dir, "inputs", "mod_inputs.Rmd")
render(pars.template.loc, output_dir = out.dir, output_format = plain_document)
render(inputs.template.loc, output_dir = out.dir, output_format = plain_document)
