## Get the variables as they were saved on the original computer
setup.env <- new.env()
Rfile <- "tmp"
load(paste0(Rfile, ".RData"), envir = setup.env)

## Want to change file locations from in.root to out.root
in.root <- "/home/phe.gov.uk/paul.birrell/Documents/PHE/stats/Wuhan_2019_Coronavirus"
## in.root <- "/project/pandemic_flu/Wuhan_Coronavirus"
## in.root <- "/rds/user/aa995/hpc-work"
if(Sys.info()["user"] == "pjb51") out.root <- "/rds/user/pjb51/hpc-work/project/pandemic_flu/Wuhan_Coronavirus"
if(Sys.info()["user"] == "jbb50") out.root <- "/home/jbb50/rds/hpc-work"

## Change location of repo
in.repo <- "real-time-mcmc"
out.repo <- "real-time-mcmc"

## Change location of output directory
in.base <- "PrevINLA424_cm6ons_IFR4bp_NHS60cutoff_18wk2_prev14-0PHE_matrices_20210709_timeuse_household_new_base_deaths"
out.base <- "NHS60_chain1"
## in.base <- "_new_base"
## out.base <- ""

## Get all variable names
var.list <- eapply(setup.env, typeof)
vbl.names <- names(var.list)[unlist(var.list) == "character"]

## Within the specified environment...
with(setup.env, {
    for(vars in vbl.names){
        assign(vars, gsub(in.root, out.root, get(vars), fixed = TRUE))
        assign(vars, gsub(in.repo, out.repo, get(vars), fixed = TRUE))
        assign(vars, gsub(in.base, out.base, get(vars), fixed = TRUE))
    }
})
## Added exception for prev.data
setup.env$prev.data <- lapply(setup.env$prev.data, function(x)
    gsub(in.root, out.root, x, fixed = TRUE))
setup.env$prev.data <- lapply(setup.env$prev.data, function(x)
    gsub(in.repo, out.repo, x, fixed = TRUE))
setup.env$prev.data <- lapply(setup.env$prev.data, function(x)
    gsub(in.base, out.base, x, fixed = TRUE))

## Temporary line to be deleted
if(exists("infections")) rm(infections)
expit <- function(x) exp(x)/(1+exp(x))
## abreaks.icr <- 3:7


save(list = ls(envir = setup.env), file = paste0(Rfile, ".RData"), envir = setup.env)
