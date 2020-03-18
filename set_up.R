require(rmarkdown)

date.of.runs <- "20200324"

## ## Regions under study
regions <- c("ENGLAND")
## Which intervention scenario
scenario.name <- "med"
## Give the run a name - will be used to save a directory for outputs.
out.dir <- paste0("initial_run_deaths_delaysensENG", date.of.runs, "_", scenario.name, "/")

source("set_up_inputs.R")
source("set_up_pars.R")


pars.template.loc <- file.path(proj.dir, "inputs", "mod_pars.Rmd")
inputs.template.loc <- file.path(proj.dir, "inputs", "mod_inputs.Rmd")
render(pars.template.loc, output_dir = out.dir, output_format = plain_document)
render(inputs.template.loc, output_dir = out.dir, output_format = plain_document)

setwd(out.dir)
system("html2text -width 9999 -o ./mod_pars.txt ./mod_pars.html")
system("html2text -width 9999 -o ./mod_inputs.txt ./mod_inputs.html")

system("nice -19 ./rtm_gnu > runtime.txt", intern = TRUE)

setwd(cur.dir)
