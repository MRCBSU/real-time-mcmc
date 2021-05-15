## Get the variables as they were saved on the original computer
setup.env <- new.env()
load("tmp.RData", envir = setup.env)

## Want to change file locations from in.root to out.root
in.root <- "/home/phe.gov.uk/paul.birrell/Documents/PHE/stats/Wuhan_2019_Coronavirus"
## in.root <- "/project/pandemic_flu/Wuhan_Coronavirus"
if(Sys.info()["user"] == "pjb51") out.root <- "/rds/user/pjb51/hpc-work/project/pandemic_flu/Wuhan_Coronavirus"
if(Sys.info()["user"] == "jbb50") out.root <- "/home/jbb50/rds/hpc-work"

## Change location of repo
in.repo <- "real-time-mcmc"
out.repo <- "real-time-mcmc-amgs"

## Change location of output directory
in.base <- "Prev354_cm4ons_IFR3bp_ONS60cutoff_25wk2_prev14-5Jamie_matrices_20210423_timeuse_household_deaths"
out.base <- file.path("Prev354_cm4ons_IFR3bp_ONS60cutoff_25wk2_prev14-5Jamie_matrices_20210423_timeuse_household_deaths", "projections_long_endpoint")

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


save(list = ls(envir = setup.env), file = "tmp.RData", envir = setup.env)
