require(tidyverse)
require(lubridate)
require(rmarkdown)
require(abind)
require(parallel)
require(knitr)

out.dir <- "./" ## commandArgs(trailingOnly = TRUE)[1]
QUANTILES <- c(0.025, 0.5, 0.975)
## out.dir <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.na(out.dir)) setwd(out.dir)
cat(getwd(), "\n")
cat(out.dir, "\n")
load("mcmc.RData")
load("tmp.RData")
source(file.path(Rfile.loc, "sim_func.R"))

projections.basename <- "projections_snapshot"
projections.basedir <- file.path(out.dir, projections.basename)

## Enter the date for which we need the snapshot
snap.date <- ymd("20211208")


## ---------------- IDEALLY CODE BELOW HERE SHOULD NOT BE CHANGED BETWEEN RUNS
if(!file.exists(projections.basedir))
    dir.create(projections.basedir)

## ## DATA STRUCTURE CHANGES TO VARIABLES REQUIRED IN THE mod_inputs.Rmd TEMPLATE
ndays <- snap.date - start.date + 1
start.hosp <- 1
start.gp <- 1
start.prev <- 1
end.hosp <- ifelse(hosp.flag, ndays, 1)
end.gp <- ifelse(gp.flag, ndays, 1)
end.prev <- ifelse(prev.flag, ndays, 1)
end.vac <- ifelse(vacc.flag, ndays, 1)

## REMOVE ANY UNNECESSARY CONTACT MATRICES
cm.flags <- c(TRUE, cm.breaks < ndays)
cm.breaks <- cm.breaks[cm.flags[-1]]
cm.bases <- cm.bases[cm.flags]
cm.mults <- cm.mults[cm.flags]

if(!all(file.exists(cm.mults)))
    stop("Specified multiplier matrix doesn't exist")

## DATA VARIABLES (and data files)
repeat.last.row <- function(real.fl, dummy.fl){
    tmpdata <- read_tsv(real.fl, col_names = FALSE)
    dummy.fl <- file.path(projections.basedir, paste0(dummy.fl, ".txt"))
    tmpdata <- bind_rows(tmpdata,
                         tmpdata[rep(nrow(tmpdata), ndays - nrow(tmpdata)), ]) %>%
            mutate(Date = start.date - 1 + (1:ndays)) %>%
            select(-X1) %>%
            select(Date, everything()) %>%
            write_tsv(dummy.fl, col_names = FALSE)
    return(dummy.fl)
}

if(gp.flag)
    for(reg in regions){
        cases.files[reg] <- repeat.last.row(cases.files[reg], paste0("dummy_cases_", reg))
        denoms.files[reg] <- repeat.last.row(denoms.files[reg], paste0("dummy_denoms_", reg))
    }
if(hosp.flag)
    for(reg in regions)
        hosp.data[reg] <- repeat.last.row(hosp.data[reg], paste0("dummy_deaths_", reg))
if(prev.flag){
    if (is.null(names(prev.data$lmeans))) {
        names(prev.data$lmeans) <- regions
        names(prev.data$lsds) <- regions
    }
    for(reg in regions){
        prev.data$lmeans[reg] <- repeat.last.row(prev.data$lmeans[reg], paste0("dummy_prev_lmeans_", reg))
        prev.data$lsds[reg] <- repeat.last.row(prev.data$lsds[reg], paste0("dummy_prev_lsds_", reg))
    }
}

## MCMC control
num.iterations <- 1
thin.params <- 1
thin.outputs <- 1
adaptive.phase <- 0
burnin <- 0
num.threads <- 1
if(projections.basename == "projections_snapshot") mcmc.outs <- smc.outs <- 1 else smc.outs <- 0

## The mod_inputs.txt file wont change with each projections so can render it now
knit(input = inputs.template.loc, output = file.path(projections.basedir, "mod_inputs.txt"))

## ## ### --- End of setting inputs for mod_inputs.txt --- ### ## ##

## ## CHANGES NECESSARY FOR THE mod_pars.Rmd TEMPLATE

## Will need to remove temporal breakpoints and rows from design matrices


## ## ### ---  End of setting inputs for mod_pars.txt  --- ### ## ##
