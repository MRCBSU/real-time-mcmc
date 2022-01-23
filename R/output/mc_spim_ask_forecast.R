require(tidyverse)
require(lubridate)
require(rmarkdown)
require(abind)
require(parallel)
require(knitr)

out.dir <- commandArgs(trailingOnly = TRUE)[1]
setwd(out.dir)
out.ind <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
Rscenario <- c(0.7, 0.9, 1.1, 1.3)[out.ind]
## Rscenario <- 0.8
QUANTILES <- c(0.025, 0.5, 0.975)

load("mcmc.RData")
load("tmp.RData")

Renv <- new.env()
load("output_matrices.RData", envir = Renv)
Rt <- Renv$Rt
## Fix the iteration numbers that we're going to look at
int_iter <- 0:(num.iterations - 1)
## parameter.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.params)
parameter.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= burnin]
## outputs.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.outputs)
outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]
parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)
stopifnot(length(parameter.to.outputs) == length(outputs.iterations)) # Needs to be subset
iterations.for.Rt <- parameter.to.outputs[seq(from = 1, to = length(parameter.to.outputs), length.out = 500)]
outputs.for.Rt <- which(parameter.to.outputs %in% iterations.for.Rt)
rm(Renv)

source(file.path(Rfile.loc, "sim_func.R"))

nweeks.ahead <- 9

counterfactual <- FALSE

projections.basedir <- file.path(out.dir, paste0("projections_MTP_R_", Rscenario))

nforecast.weeks <- nweeks.ahead - nforecast.weeks
mm.breaks <- start.date - 1 + max(cm.breaks) + (1:nforecast.weeks * days(7))
google.data.date <- ymd(google.data.date)
mult.order <- rep(1, length(mm.breaks))
sero.flag <- 0 ## Are we interested in simulating serological outputs? Switched off for the moment.
prev.flag <- 1 ## prev.flag ## Are we interested in simulating prevalence outputs?

if(prev.flag && any(prev.data$lmeans == "NULL")){
    ## if (!exists("date.prev")) {
		## Get the date of the prevalence data
		date.prev <- ymd("20210317")
		## Convert that to an analysis day number
		prev.end.day <- 388
		last.prev.day <- 388
		first.prev.day <- 75
		if(!exists("days.between.prev")) days.between.prev <- 7
		## Default system for getting the days on which the likelihood will be calculated.
		prev.lik.days <- rev(seq(from = last.prev.day, to = first.prev.day, by = -days.between.prev))
	## }
    for(r in 1:nr){
	  prev.file.prefix <- paste0(data.dirs["prev"], "/", date.prev, "_", paste(prev.lik.days, collapse = "_"), "_", regions[r], "ons_") ## , paste0(prev.lik.days, collapse = "_"), "_")
          prev.file.suffix <- paste0("logprev.txt")
      prev.data$lmeans[r] <- paste0(prev.file.prefix, "mean", prev.file.suffix)
      prev.data$lsds[r] <- paste0(prev.file.prefix, "sd", prev.file.suffix)
    }
    names(prev.data$lmeans) <- names(prev.data$lsds) <- regions
}
## ## Do the contract matrices that will be used in the projection exist on file and are they correct?
overwrite.matrices <- FALSE

## ## mod_pars.Rmd specifications that will change - should only be breakpoints and design matrices
if(prev.flag & all(prior.r1 == 1)) value.r1 <- 7.18
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

ndays <- lubridate::as_date(date.data) - start.date + (7 * nweeks.ahead) + 1
start.hosp <- 1
start.gp <- 1
start.prev <- 1
end.hosp <- ifelse(hosp.flag | adm.flag, ndays, 1)
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
if(hosp.flag | adm.flag)
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
smc.outs <- 0

## The mod_inputs.txt file wont change with each projections so can render it now
## render(inputs.template.loc, output_dir = file.path(out.dir, "projections"), output_format = "plain_document")
knit(input = inputs.template.loc, output = file.path(projections.basedir, "mod_inputs.txt"))

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
if(!single.ifr)
    symlink.design("ifr.design.txt")
if(vacc.flag){
    symlink.design("vac.alpha1.design.txt")
    symlink.design("vac.alphan.design.txt")
    symlink.design("vac.pi1.design.txt")
    symlink.design("vac.pin.design.txt")
}
## The matrix for the random-walks will need changing to account for an extra breakpoint
## Place the new break-point now
today.break <- ymd("20220126") - start.date + 1
beta.breaks <- c(beta.breaks, today.break)
beta.block <- beta.design[1:nbetas, 1:nbetas] %>%
    rbind(rep(1, nbetas)) %>%
    cbind(c(rep(0, nbetas), 1))
require(Matrix)
beta.design <- lapply(1:nr, function(x) beta.block) %>%
    bdiag() %>%
    as.matrix()
write_tsv(as.data.frame(beta.design), file.path(projections.basedir, "beta.design.txt"), col_names = FALSE)
## ## ## --------------------------------------------------------------

## Subset the MCMC chain to look at only the iterations for which we have calculated an R value
params <- lapply(params, function(x) x[iterations.for.Rt, , drop = FALSE])
## Append a new beta to hit the target R value
endRt <- Rt[, min(c(dim(Rt)[2], today.break)), ]
beta.new <- log(Rscenario / endRt)

## ## ## MAIN PROJECTION LOOP

## ## Set-up output quantities
NNI <- NNI.files <- vector("list", nr)
Sero <- Sero.files <- vector("list", nr)
if(vacc.flag) DNNI <- DNNI.files <- vector("list", nr)
if(hosp.flag) Deaths <- Deaths.files <- vector("list", nr)
if(gp.flag) Cases <- Cases.files <- vector("list", nr)
if(prev.flag) Prevs <- Prev.files <- vector("list", nr)
params$log_beta_rw <- params$log_beta_rw %>%
    array(dim = c(dim(params$log_beta_rw)[1], nbetas, nr)) %>%
    abind(beta.new, along = 2) %>%
    apply(1, as.vector) %>%
    t()
nbetas <- nbetas + 1

## Get iteration indices
niter <- length(iterations.for.Rt)

## ## For each iteration
pct <- 0
## xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1)
if(Sys.info()["user"] %in% c("jbb50", "pjb51")){
    exe <- "hpc2"
} else exe <- Sys.info()["nodename"]
cat("rtm.exe = ", exe, "\n")
cat("full file path = ", file.path(proj.dir, paste0("../real-time-mcmc-dev/rtm_", exe)), "\n")
proj.dir <- gsub("amgs", "dev", proj.dir, fixed = TRUE)
xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1, rtm.exe = exe)
proj.dir <- gsub("dev", "amgs", proj.dir, fixed = TRUE)
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
dim.list <- list(iteration = iterations.for.Rt,
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
if(hosp.flag | adm.flag) {
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
save(list = save.list, file = paste0("projections_R", Rscenario, ".RData"))

## ## ## Housekeeping
lapply(hosp.data, file.remove)
lapply(cases.files, file.remove)
lapply(denoms.files, file.remove)
if(prev.flag)
    lapply(prev.data, file.remove)
