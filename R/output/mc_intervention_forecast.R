require(tidyverse)
require(lubridate)
require(rmarkdown)
require(abind)
require(parallel)
require(knitr)

load("mcmc.RData")
load("tmp.RData")

source(file.path(Rfile.loc, "sim_func.R"))

QUANTILES <- c(0.025, 0.5, 0.975)
array.num <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## ## mod_inputs.Rmd items that will change in the projections.

## Forecast projection
nforecast.weeks <- 24

intervention.scenario <- ifelse(array.num <= 4, "realistic", "optimistic")
if (array.num %% 2 == 0) intervention.scenario <- paste0("school_", intervention.scenario)
intervention.length.weeks <- ifelse(array.num %in% c(1, 2, 5, 6), 2, 4)
intervention.start <- ymd(20201105)
intervention.matrix <- paste0(intervention.scenario, "2020-11-05.csv")
projections.file <- paste0("projections_", intervention.scenario, "_", intervention.length.weeks, "weeks.RData")
projections.basedir <- file.path(out.dir, paste0("projections", array.num))

## Enter dates at which it is anticipated that the contact model will change
## mm.breaks <- ymd("20200928") + (1:nforecast.weeks * days(7))
mm.breaks <- cm.breaks[length(cm.breaks)] + (1:nforecast.weeks) * 7 + start.date - 1
google.data.date <- ymd("20201030")
intervention.breaks <- c(intervention.start, intervention.start + intervention.length.weeks * 7)
mult.order <- rep(1, length(mm.breaks))

## ## ----------------------------------------------------------

## ## mod_pars.Rmd specifications that will change - should only be breakpoints and design matrices

bank.holiday.days.new <- NULL
## ## ---------------------------------------------------------------------------------------------

## ## ## FUNCTION DEFINITIONS
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
symlink.design <- function(design)
    file.symlink(file.path(out.dir, design), projections.basedir)
## ## Compile a forecast output
combine.rtm.output <- function(x, strFld){
    oList <- lapply(x, function(x) do.call(abind, args = list(x[[strFld]], along = 3)))
    oList <- do.call(abind, args = list(oList, along = 0))

    }

## ## ## --------------------

if(!file.exists(projections.basedir))
    dir.create(projections.basedir)

## ## ## CHANGES TO VARIABLES BASED ON mod_inputs-LIKE SPECIFICATIONS
sero.flag <- 0 ## Are we interested in serological outputs? Switched off for the moment.
ndays <- lubridate::as_date(date.data) - start.date + (7 * nforecast.weeks) + 1
start.hosp <- 1
start.gp <- 1
end.hosp <- ifelse(hosp.flag, ndays, 1)
end.gp <- ifelse(gp.flag, ndays, 1)
prev.flag <- 0

## Get the new contact matrices to use
cm.breaks <- c(cm.breaks, mm.breaks - start.date + 1)
intervention.breaks <- intervention.breaks - start.date + 1
idx <- which(cm.breaks %in% seq(intervention.breaks[1], intervention.breaks[2]))
cm.files <- c(cm.files,
              paste0("england_8ag_contact_ltprojwk", 1:length(mm.breaks), "_", format(google.data.date, "%Y%m%d"), ".txt"))
cm.bases <- file.path(proj.dir, "contact_mats", cm.files)
## Use the next line to specify where the new matrices are stored
matrix.dir <- file.path(dirname(matrix.dir), paste0("google_mobility_relative_matrices_", format(google.data.date, "%Y%m%d")))
intervention.dir <- file.path(dirname(matrix.dir), paste0("scenarios_", format(google.data.date, "%Y%m%d")))
## Use the next line to specify the format with the filenames
## cm.lockdown.fl <- c(cm.lockdown.fl, paste0("England", mm.breaks, "all.csv"))
cm.lockdown.fl <- c(cm.lockdown.fl, paste0("England", mm.breaks, "all.csv"))
cm.intervention.fl <- file.path(intervention.dir, intervention.matrix)
cm.lockdown <- c(cm.lockdown,
                 file.path(matrix.dir, tail(cm.lockdown.fl, length(mm.breaks))))

cm.breaks <- c(cm.breaks[1:(first(idx) - 1)], intervention.breaks, cm.breaks[(last(idx)+1):length(cm.breaks)])
cm.lockdown <- c(cm.lockdown[1:(first(idx) - 1)], cm.intervention.fl, cm.lockdown[last(idx):length(cm.lockdown)])
cm.bases <- c(cm.bases[1:first(idx)],
	    file.path(proj.dir, "contact_mats", paste0("england_8ag_contact_intervention", "_", intervention.scenario, format(google.data.date, "%Y%m%d"), ".txt")),
            cm.bases[(last(idx) + 1):length(cm.bases)])

if(!all(file.exists(cm.bases))){
    idx.miss <- which(!file.exists(cm.bases))
    adf <- as.data.frame(lst$England$all$m)
    for(idx in idx.miss){
        mat <- (read_csv(cm.lockdown[idx - 1]) * adf) %>%
            write_tsv(cm.bases[idx], col_names = FALSE)
    }
}
stopifnot(all(file.exists(cm.bases)))
## Get the new contract matrix modifiers to use
cm.mults <- c(cm.mults,
              file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_mod3levels", mult.order, ".txt"))
              )[1:length(cm.bases)]
if(!all(file.exists(cm.mults)))
    stop("Specified multiplier matrix doesn't exist")
stopifnot(length(cm.bases) == length(cm.mults))
stopifnot(length(cm.bases) == length(cm.breaks) + 1)

## Need some dummy data files padded to the right number of rows
if(gp.flag)
    for(reg in regions){
        cases.files[reg] <- repeat.last.row(cases.files[reg], paste0("dummy_cases_", reg))
        denoms.files[reg] <- repeat.last.row(denoms.files[reg], paste0("dummy_denoms_", reg))
    }
if(hosp.flag)
    for(reg in regions)
        hosp.data[reg] <- repeat.last.row(hosp.data[reg], paste0("dummy_deaths_", reg))
if(prev.flag)
    for(reg in regions){
        prev.data$lmeans[reg] <- repeat.last.row(prev.data$lmeans[reg], paste0("dummy_prev_lmeans_", reg))
        prev.data$lsds[reg] <- repeat.last.row(prev.data$lsds[reg], paste0("dummy_prev_lsds_", reg))
    }

if(!exists("study_region_str")) study_region_str <- ""

## MCMC control
num.iterations <- 1
thin.params <- 1
thin.outputs <- 1
adaptive.phase <- 0
burnin <- 0
num.threads <- 1

## The mod_inputs.txt file wont change with each projections so can render it now
## render(inputs.template.loc, output_dir = file.path(out.dir, "projections"), output_format = plain_document())
knit(input = inputs.template.loc, output = file.path(projections.basedir, "mod_inputs.txt"))

## ## ## ------------------------------------------------------------

## ## ## CHANGES TO VARIABLES BASED ON mod_pars.txt-LIKE SPECIFICATIONS

## Extend the breakpoints for the day of week effects
if(gp.flag){
    require(Matrix)
    bank.holiday.days <- c(bank.holiday.days, bank.holiday.days.new)
    ll.days <- start.date + days(0:(ndays-1))
    DAYS <- wday(ll.days, label = TRUE)
    DAYS[ll.days %in% bank.holiday.days] <- "Sun"
    lm.mat <- model.matrix(~DAYS, contrasts = list(DAYS = "contr.sum"))[, -1] %>%
        as.data.frame() %>%
        write_tsv(file.path(projections.basedir, "d_o_w_design_file.txt"), col_names = FALSE)
}

## Copy to projection directory other design matrices
if(gp.flag)
    symlink.design("icr.design.txt")
if(rw.flag)
    symlink.design("m.design.txt")
if(beta.rw.flag)
    symlink.design("beta.design.txt")
if(!single.ifr)
    symlink.design("ifr.design.txt")
## ## ## --------------------------------------------------------------

## ## ## MAIN PROJECTION LOOP

## ## Set-up output quantities
NNI <- NNI.files <- vector("list", nr)
if(hosp.flag) Deaths <- Deaths.files <- vector("list", nr)
if(gp.flag) Cases <- Cases.files <- vector("list", nr)
if(prev.flag) Prevs <- Prev.files <- vector("list", nr)

## ## Get number of iterations
niter <- min(sapply(params, nrow))

## ## For each iteration
pct <- 0
## xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1)
if(Sys.info()["user"] %in% c("pjb51", "jbb50")){
    exe <- "hpc"
} else exe <- Sys.info()["nodename"]
cat("rtm.exe = ", exe, "\n")
cat("full file path = ", file.path(proj.dir, paste0("rtm_", exe)), "\n")
xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1, rtm.exe = exe)

NNI <- lapply(xtmp, function(x) x$NNI)
Deaths <- lapply(xtmp, function(x) x$Deaths)
Cases <- lapply(xtmp, function(x) x$Cases)
Prevs <- lapply(xtmp, function(x) x$Prevs)
rm(xtmp)

## names(NNI) <- regions
## if(dths.flag) names(Deaths) <- regions
## if(cases.flag) names(Cases) <- regions
## ## ## --------------------

melt.list <- function(xlist)
    abind(lapply(xlist,
                 abind,
                 along = 3),
          along = 0)

## ## ## SAVE SOME OUTPUTS
dim.list <- list(iteration = 1:niter,
                 age = age.labs,
                 date = start.date + 0:(ndays - 1),
                 region = regions
                 )
infections <- melt.list(NNI);rm(NNI)
save.list <- "infections"
dimnames(infections) <- dim.list
if(hosp.flag) {
    deaths <- melt.list(Deaths);rm(Deaths)
    save.list <- c(save.list, "deaths")
    dimnames(deaths) <- dim.list
    }
if(gp.flag){
    cases <- melt.list(Cases);rm(Cases)
    save.list <- c(save.list, "cases")
    dimnames(cases) <- dim.list
}
if(prev.flag){
    prevalence <- melt.list(Prevs);rm(Prevs)
    save.list <- c(save.list, "prevalence")
    dimnames(prevalence) <- dim.list
}
save(list = save.list, file = projections.file)

## ## ## ## Housekeeping
lapply(hosp.data, file.remove)
lapply(cases.files, file.remove)
lapply(denoms.files, file.remove)
if(prev.flag)
    lapply(prev.data, file.remove)
rm(list=ls())
