require(rmarkdown)

## Give the run a name - will be used to save a directory for outputs.
out.dir <- "initial_run_linelist_egr/"

source("set_up_inputs.R")
source("set_up_pars.R")

plain_document <- output_format(
    knitr = knitr_options(),
    pandoc = pandoc_options(to = "plain", ext = ".txt"),
)

render("inputs/mod_pars.Rmd", output_dir = out.dir, output_format = plain_document)
render("inputs/mod_inputs.Rmd", output_dir = out.dir, output_format = plain_document)
setwd(out.dir)

system("./rtm_gnu > runtime.txt", intern = TRUE)
