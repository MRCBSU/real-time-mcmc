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

projections.basename <- "projections_snapshot"
projections.basedir <- file.path(out.dir, projections.basename)

## Enter the date for which we need the snapshot
snap.date <- ymd("20211026")

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
if(start.sero > ndays) sero.flag <- 0
end.sero <- ifelse(sero.flag, ndays, 1)

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
    if(ndays > nrow(tmpdata)){
    tmpdata <- bind_rows(tmpdata,
                         tmpdata[rep(nrow(tmpdata), ndays - nrow(tmpdata)), ]) %>%
            mutate(Date = start.date - 1 + (1:ndays)) %>%
            select(-X1) %>%
            select(Date, everything()) %>%
            write_tsv(dummy.fl, col_names = FALSE)
    } else write_tsv(tmpdata[1:ndays, ], dummy.fl, col_names = FALSE)
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

## ## ## ### --- End of setting inputs for mod_inputs.txt --- ### ## ## ##

## ## ## CHANGES NECESSARY FOR THE mod_pars.Rmd TEMPLATE

## ## Will need to remove temporal breakpoints and rows from design matrices
vac.rows.remove <- function(t0, nr.levels, na.levels, t.levels){
    nt.levels <- length(t.levels)
    tot.rows <- nr.levels * na.levels * nt.levels
    as_tibble(expand.grid(1:na.levels, t.levels, 1:nr.levels)) %>%
        ## which levels do we want to keep...
        mutate(include = Var2 <= t.levels[t0]) %>%
        pull(include)
}
copy.design <- function(design)
    file.copy(file.path(out.dir, design), projections.basedir)
symlink.design <- function(design)
    file.symlink(file.path(out.dir, design), projections.basedir)
## Select only relevant rows from the original design matrix according to the desired timepoint.
select.design <- function(design, breaks, t0, nreg = nr, nA = nA){
    ## read in the original design matrix
    des <- read_tsv(file.path(out.dir, design), col_names = FALSE)
    ## how many of the temporal breakpoints are actually active
    nbg <- sum(breaks < t0)
    ## which rows of the original design matrix will we now have to remove
    bg.include <- vac.rows.remove(nbg + 1, nreg, nA, 1:(length(breaks) + 1))
    ## return only the rows we want from the original design
    des[bg.include, ]
}

## ## Copy across the original design matrices
if(gp.flag)
    symlink.design("icr.design.txt")
if(rw.flag)
    symlink.design("m.design.txt")
if(beta.rw.flag){
    ## Can we use the pre-existing matrix?
    if(max(beta.breaks) < ndays){ ## yes
        symlink.design("beta.design.txt")
    } else {
        select.design("beta.design.txt", beta.breaks, ndays, nA = 1) %>%
            write_tsv(path = file.path(projections.basedir, "beta.design.txt"), col_names = FALSE)
        beta.breaks <- beta.breaks[beta.breaks < ndays]
        ## des <- read_tsv(file.path(out.dir, "beta.design.txt"), col_names = FALSE)
        ## nbg <- sum(beta.breaks < ndays)
        ## bg.include <- vac.rows.remove(nbg + 1, nr, 1, 1:(length(beta.breaks) + 1))
        ## write_tsv(des[bg.include, ], path = file.path(projections.basedir, "beta.design.txt"), col_names = FALSE)
    }
}
if(!single.ifr){
    ## Can we use the pre-existing matrix?
    if(max(tbreaks.ifr) < ndays){
        symlink.design("ifr.design.txt")
    } else {
        select.design("ifr.design.txt", tbreaks.ifr, ndays, nr = 1, nA = nA - 1) %>%
            write_tsv(path = file.path(projections.basedir, "ifr.design.txt"), col_names = FALSE)
        tbreaks.ifr <- tbreaks.ifr[tbreaks.ifr < ndays]
    }
}
if(vacc.flag){
    ## Can we use the pre-existing design matrix?
    if(max(vac.effec.bp) < ndays){
        symlink.design("vac.pi1.design.txt")
        symlink.design("vac.pin.design.txt")
        symlink.design("vac.alpha1.design.txt")
        symlink.design("vac.alphan.design.txt")
    } else {
        mat1 <- select.design("vac.pi1.design.txt", vac.effec.bp, ndays, nr, nA)
        matn <- select.design("vac.pin.design.txt", vac.effec.bp, ndays, nr, nA)
        write_tsv(mat1, path = file.path(projections.basedir, "vac.pi1.design.txt"), col_names = FALSE)
        write_tsv(mat1, path = file.path(projections.basedir, "vac.alpha1.design.txt"), col_names = FALSE)
        write_tsv(matn, path = file.path(projections.basedir, "vac.pin.design.txt"), col_names = FALSE)
        write_tsv(matn, path = file.path(projections.basedir, "vac.alphan.design.txt"), col_names = FALSE)
        vac.effec.bp <- vac.effec.bp[vac.effec.bp < ndays]
    }
}
## ## -------------------------------------

## ## ## ### ---  End of setting inputs for mod_pars.txt  --- ### ## ## ##

## ## ## ### --- --- ---   MAIN PROJECTION LOOP   --- --- --- ### ## ## ##

## ## Get number of iterations
niter <- min(sapply(params, nrow))

## ## For each iteration
pct <- 0
## xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1)
if(Sys.info()["user"] %in% c("jbb50", "pjb51")){
    exe <- "hpc2"
} else exe <- Sys.info()["nodename"]
cat("rtm.exe = ", exe, "\n")
xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 1, rtm.exe = exe)
NNI <- lapply(xtmp, function(x) x$NNI)
if(smc.outs){
    state <- lapply(xtmp, function(x) x$state)
} else {
    Sero <- lapply(xtmp, function(x) x$Sero)
    if(vacc.flag) DNNI <- lapply(xtmp, function(x) x$DNNI)
    Deaths <- lapply(xtmp, function(x) x$Deaths)
    Cases <- lapply(xtmp, function(x) x$Cases)
    Prevs <- lapply(xtmp, function(x) x$Prevs)
}
rm(xtmp)
## ## ## ###   --- --- --- --- --- --- --- --- --- --- ---   ### ## ## ##


## ## PROCESSING OUTPUTS
melt.list <- function(xlist)
    abind(lapply(xlist,
                 abind,
                 along = 3),
          along = 0)
dim.list <- list(iteration = 1:niter,#niter,
                 age = age.labs,
                 date = paste0("t", seq(0, as.integer(ndays) - 0.5, by = 0.5)),
                 region = regions
                 )

infections <- melt.list(NNI)## ;rm(NNI)
dimnames(infections) <- dim.list
state <- do.call(bind_rows, state) %>%
    inner_join(expand_grid(regions, age.labs) %>% mutate(popn = pop.input) %>% rename(region = regions))
state <- bind_rows(state,
                   state %>% group_by(iteration, region, age.labs) %>%
                   summarise(unrecovered = sum(value), popn = median(popn)) %>%
                   mutate(value = popn - unrecovered,
                          state.names = "R-") %>%
                   select(-unrecovered)
                   )
state.lkup <- tibble(state.names = c("S", "SV1", "SV2", "E1", "E2", "I1", "I2", "R+", "R-", "plambda"), state.text = c("Fully susceptible", "Never infected, incomplete immunisation", "Never infected, complete immunisation", "Latent infection", "Latent infection", "Prevalent infection", "Prevalent infection", "Prevalent infection", "Infection-acquired immunity", "Infectious pressure"))
state <- state %>%
    inner_join(state.lkup)

save(state, file = paste0("snapshot", snap.date, ".RData"))
