require(tidyverse)
require(lubridate)
require(rmarkdown)
require(abind)
require(parallel)
require(knitr)

out.dir <- commandArgs(trailingOnly = TRUE)[1]
QUANTILES <- c(0.025, 0.5, 0.975)
## out.dir <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.na(out.dir)) setwd(out.dir)
cat(getwd(), "\n")
cat(out.dir, "\n")
load("mcmc.RData")
load("tmp.RData")
source(file.path(Rfile.loc, "sim_func.R"))
## ## mod_inputs.Rmd items that will change in the projections.
## Number of weeks to forecast ahead
nweeks.ahead <-11

counterfactual <- FALSE



projections.basename <-"projections_counter"
#projections.basename <- "projections_midterm"




projections.basedir <- file.path(out.dir, projections.basename)
## ## Enter dates at which it is anticipated that the contact model will change
## mm.breaks <- ymd("20201109") + (1:nforecast.weeks * days(7))
## ## Forecast projection
nforecast.weeks <- nweeks.ahead - nforecast.weeks
mm.breaks <- start.date - 1 + max(cm.breaks) + (1:nforecast.weeks * days(7))
google.data.date <- ymd(google.data.date)
mult.order <- rep(1, length(mm.breaks))
sero.flag <- 0 ## Are we interested in simulating serological outputs? Switched off for the moment.
if(prev.flag && any(prev.data$lmeans == "NULL")){
	stop("No prev files specified")
}
prev.flag <- 1 ## prev.flag ## Are we interested in simulating prevalence outputs?
if(prev.flag && any(prev.data$lmeans == "NULL")){
#This HEAD up to COVID_vacc_amgs was the suggestion from Colin's/Paul's branch but seems that my/Josh's version is proper to use.
#<<<<<<< HEAD
#=======
    ## if (!exists("date.prev")) {
		## Get the date of the prevalence data








		date.prev <- ymd("20210324")
		## Convert that to an analysis day number
		prev.end.day <- 395
		last.prev.day <- 395
		first.prev.day <- 75
		if(!exists("days.between.prev")) days.between.prev <- 7

		## Default system for getting the days on which the likelihood will be calculated.
#		prev.lik.days <- rev(seq(from = last.prev.day, to = first.prev.day, by = -days.between.prev))
	## }
#>>>>>>> COVID_vacc_amgs
    for(r in 1:nr){
      prev.data$lmeans[r] <- file.path(data.dirs["prev"], "dummy_meanlogprev.txt")
      prev.data$lsds[r] <- file.path(data.dirs["prev"], "dummy_sdlogprev.txt")
    }
    names(prev.data$lmeans) <- names(prev.data$lsds) <- regions
}
## ## Do the contract matrices that will be used in the projection exist on file and are they correct?
overwrite.matrices <- FALSE

## ## ----------------------------------------------------------

## ## mod_pars.Rmd specifications that will change - should only be breakpoints and design matrices
if(prev.flag & all(prior.r1 == 1)) value.r1 <- 7.18
if(grepl("projections_counter", projections.basedir, fixed = TRUE))
{
    value.vac.alpha1 <- rep(0, length(value.vac.alpha1))
    value.vac.alphan <- rep(0, length(value.vac.alpha2))
    value.vac.pi1 <- rep(0, length(value.vac.pi1))
    value.vac.pi2 <- rep(0, length(value.vac.pi2))
}
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
ndays <- lubridate::as_date(date.data) - start.date + (7 * nweeks.ahead) + 1
start.hosp <- 1
start.gp <- 1
start.prev <- 1
end.hosp <- ifelse(hosp.flag, ndays, 1)
end.gp <- ifelse(gp.flag, ndays, 1)
end.prev <- ifelse(prev.flag, ndays, 1)
end.vac <- ifelse(vacc.flag, ndays, 1)

## Get the new contact matrices to use
cm.breaks <- c(cm.breaks, mm.breaks - start.date + 1)
cm.files <- c(cm.files,
              paste0("england_8ag_contact_projwk", 1:length(mm.breaks), "_", google.data.date.str, ".txt"))
cm.bases <- file.path(proj.dir, "contact_mats", cm.files)
cm.lockdown.fl <- c(cm.lockdown.fl, paste0("England", mm.breaks, "all.csv"))
cm.lockdown <- c(cm.lockdown,
                 file.path(matrix.dir, tail(cm.lockdown.fl, length(mm.breaks))))
## Get the new contract matrix modifiers to use
cm.mults <- c(cm.mults,
              file.path(proj.dir,
                        "contact_mats",
                        paste0("ag", nA, "_mult_mod", ifelse(contact.model!=6, contact.model, "All"), "Levels", mult.order, ".txt"))
              )
if(counterfactual){
    cm.dates <- start.date + cm.breaks - 1
    future.mats <- cm.dates > google.data.date
    cm.breaks <- cm.breaks[!future.mats]
    cm.bases <- cm.bases[c(TRUE, !future.mats)]
    cm.mults <- cm.mults[c(TRUE, !future.mats)]
}
## Write to file any contact matrix that doesn't already exist.
if(!all(file.exists(cm.bases)) || overwrite.matrices){
    idx.miss <- which(!file.exists(cm.bases))
    if(overwrite.matrices) idx.miss <- 2:length(cm.bases)
    adf <- as.data.frame(lst$England$all$m)
    for(idx in idx.miss){
        mat <- (read_csv(cm.lockdown[idx - 1]) * adf) %>%
            write_tsv(cm.bases[idx], col_names = FALSE)
    }
}
if(!all(file.exists(cm.mults)))
    stop("Specified multiplier matrix doesn't exist")

## Need some dummy data files padded to the right number of rows
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

## The mod_inputs.txt file wont change with each projections so can render it now
## render(inputs.template.loc, output_dir = file.path(out.dir, "projections"), output_format = "plain_document")
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
if(vacc.flag){
    symlink.design("vac.pi1.design.txt")
    symlink.design("vac.pin.design.txt")
    symlink.design("vac.alpha1.design.txt")
    symlink.design("vac.alphan.design.txt")
    symlink.design("vac.pi1.design.txt")
    symlink.design("vac.pin.design.txt")
}
## ## ## --------------------------------------------------------------

## ## ## MAIN PROJECTION LOOP

## ## Set-up output quantities
NNI <- NNI.files <- vector("list", nr)
Sero <- Sero.files <- vector("list", nr)
if(vacc.flag) DNNI <- DNNI.files <- vector("list", nr)
if(hosp.flag) Deaths <- Deaths.files <- vector("list", nr)
if(gp.flag) Cases <- Cases.files <- vector("list", nr)
if(prev.flag) Prevs <- Prev.files <- vector("list", nr)

## ## Get number of iterations
niter <- min(sapply(params, nrow))

## ## For each iteration
pct <- 0
## xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1)

if(Sys.info()["user"] %in% c("jbb50", "pjb51","aa995")){
    exe <- "hpc2"
} else exe <- Sys.info()["nodename"]
cat("rtm.exe = ", exe, "\n")
cat("full file path = ", file.path(proj.dir, paste0("../real-time-mcmc-amgs/rtm_", exe)), "\n")
xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1, rtm.exe = exe)
NNI <- lapply(xtmp, function(x) x$NNI)
Sero <- lapply(xtmp, function(x) x$Sero)
if(vacc.flag) DNNI <- lapply(xtmp, function(x) x$DNNI)
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
dimnames(infections) <- dim.list
seropos <- melt.list(Sero);rm(Sero)
dimnames(seropos) <- dim.list
save.list <- c("infections", "seropos")
if(vacc.flag) {
    vacc.infections <- melt.list(DNNI);rm(DNNI)
    save.list <- c(save.list, "vacc.infections")
    dimnames(vacc.infections) <- dim.list
}
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
save(list = save.list, file = paste0(projections.basename, ".RData"))

## ## ## Housekeeping
lapply(hosp.data, file.remove)
lapply(cases.files, file.remove)
lapply(denoms.files, file.remove)
if(prev.flag)
    lapply(prev.data, file.remove)
