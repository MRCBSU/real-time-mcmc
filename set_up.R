require(rmarkdown)

## Give the run a name - will be used to save a directory for outputs.
out.dir <- "initial_run_deaths_LondonSep20200319/"

source("set_up_inputs.R")
source("set_up_pars.R")

render("inputs/mod_pars.Rmd", output_dir = out.dir)
render("inputs/mod_inputs.Rmd", output_dir = out.dir)
setwd(out.dir)
system("html2text -width 9999 -o ./mod_pars.txt ./mod_pars.html")
system("html2text -width 9999 -o ./mod_inputs.txt ./mod_inputs.html")

system("nice -19 ./rtm_gnu > runtime.txt", intern = TRUE)

setwd(cur.dir)
