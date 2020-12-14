require(tidyverse)
require(lubridate)
require(rmarkdown)
require(abind)
require(parallel)
require(knitr)

out.dir <- commandArgs(trailingOnly = TRUE)[1]
## out.dir <- getwd()
QUANTILES <- c(0.025, 0.5, 0.975)

setwd(out.dir)
load("mcmc.RData")
load("tmp.RData")
source(file.path(Rfile.loc, "sim_func.R"))
array.num <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
## ## mod_inputs.Rmd items that will change in the projections.
## Forecast projection
nweeks.ahead <- 9
nforecast.weeks <- nweeks.ahead - nforecast.weeks
ndays <- lubridate::as_date(date.data) - start.date + (7 * nweeks.ahead) + 1

## save.image("tempproj.RData")
## stop()

google.data.date <- ymd(google.data.date)
cm.file.base <- "england_8ag_contact"
scenario.names <- c("low-low", "medium-low", "high-low", "low-high", "medium-high", "high-high")
scenario.name <- scenario.names[array.num]## paste0(scenario.transmission, "-", scenario.christmas)
scenario.num <- array.num## which(scenario.names == scenario.name)
scenario.transmission <- unlist(str_split(scenario.name, "-"))[1]
scenario.christmas <- unlist(str_split(scenario.name, "-"))[2]
prev.flag <- 0
if(prev.flag & (prev.data$lmeans == "NULL")){
    for(r in 1:nr){
        prev.data$lmeans[r] <- file.path(data.dirs["prev"], paste0("2020-12-02_", regions[r], "_ons_meanlogprev_286every14.txt"))
        prev.data$lsds[r] <- file.path(data.dirs["prev"], paste0("2020-12-02_", regions[r], "_ons_sdlogprev_286every14.txt"))
    }
    names(prev.data$lmeans) <- names(prev.data$lsds) <- regions
}

intervention.matrix <- function(name) {
  ## Use the next line to specify where the new matrices are stored
  intervention.dir <- file.path(dirname(matrix.dir), paste0("scenarios_", format(google.data.date, format = "%Y%m%d")))
  return(file.path(intervention.dir, name))
}

intervention.matrix.out <- function(name) {
  file.path(proj.dir, "contact_mats", paste0(cm.file.base, "_", name, "_", format(google.data.date, format = "%Y%m%d"), ".txt"))
}

use.previous.matrix <- function(new.date, old.date) {
  stopifnot(old.date %in% names(matrix.sources))
  matrix.sources[[new.date]] <<- matrix.sources[[old.date]]
  matrix.used[[new.date]] <<- matrix.used[[old.date]]
}

remove.matrices.after <- function(date) {
  matrix.sources <<- matrix.sources[ymd(names(matrix.sources)) <= date]
  matrix.used <<- matrix.used[ymd(names(matrix.used)) <= date]
 }

matrix.series <- function(last.date.used, after.which.date, last.date, start.week){
    dates <- seq(last.date.used, by = 7, to = last.date)
    weeks <- 0:(length(dates) - 1)
    inds <- dates >= after.which.date
    dates <- dates[inds]
    weeks <- weeks[inds]
    for(dt in dates){
        dta <- lubridate::as_date(dt)
        wka <- weeks[dates == dt]
        matrix.sources[[as.character(dta)]] <<- gsub(last.date.used, dta, matrix.sources[[as.character(last.date.used)]], fixed = TRUE)
        matrix.used[[as.character(dta)]] <<- gsub(paste0("ldwk", start.week), paste0("projwk", wka), matrix.used[[as.character(last.date.used)]], fixed = TRUE)
    }
}

matrix.sources <- as.list(cm.lockdown)
names(matrix.sources) <- c(cm.breaks) + start.date - 1
matrix.used <- as.list(cm.bases)
names(matrix.used) <- c(1, cm.breaks) + start.date - 1
last.forecast.date <- as.Date(last(names(matrix.used)), format = "%Y-%m-%d")
last.forecast.week <- as.integer(sub(".*?ldwk.*?(\\d+).*", "\\1", last(matrix.used)))

if(scenario.transmission != "low") lockdown.end <- ymd(20201202) else lockdown.end <- NULL; str.lockdown.end <- as.character(lockdown.end)
if(scenario.transmission != "low") holidays.end <- ymd(20210104) else holidays.end <- NULL; str.holidays.end <- as.character(holidays.end)

holidays.start <- ymd(20201221); str.holidays.start <- as.character(holidays.start)
tiers.end <- ymd(20201222); str.tiers.end <- as.character(tiers.end)
tiers.begin <- ymd(20201227); str.tiers.begin <- as.character(tiers.begin)
## periodic.breaks <- rev(seq(holidays.start, lockdown.end, by = -14)); str.periodic.breaks <- as.character(periodic.breaks)
intervention.breaks <- c(tiers.end, tiers.begin)
if(!is.null(lockdown.end)) intervention.breaks <- c(lockdown.end, intervention.breaks)
if(!is.null(holidays.end)) intervention.breaks <- c(intervention.breaks, holidays.end)

## Remove projection matrices to be over-written by the new period
remove.matrices.after(min(intervention.breaks))

## Add in the new matrix inputs
if(!is.null(lockdown.end)){
    desc.str <- ifelse(scenario.transmission == "medium", "medtrans", "hightrans")
    matrix.sources[[str.lockdown.end]] <- intervention.matrix(paste0(desc.str, intervention.breaks[1], ".csv"))
    ## Add in the new matrix outputs
    matrix.used[[str.lockdown.end]] <- intervention.matrix.out(paste0(desc.str, "_mat1"))
    which.mat <- which(grepl("2020-12-21", cm.lockdown))
    matrix.sources[[str.holidays.start]] <- cm.lockdown[which.mat]
    matrix.used[[str.holidays.start]] <- cm.bases[which.mat + 1]
}
## for(strDate in str.periodic.breaks){
##     matrix.sources[[strDate]] <- intervention.matrix(paste0("current", strDate, ".csv"))
##     matrix.used[[strDate]] <- intervention.matrix.out(paste0("current_mat", which(str.periodic.breaks %in% strDate)))
## }
if(scenario.christmas == "low"){
    matrix.sources[[str.tiers.end]] <- intervention.matrix(paste0("current", last.forecast.date, ".csv"))
    matrix.used[[str.tiers.end]] <- intervention.matrix.out("current_mat")
} else if(scenario.christmas == "high"){
    matrix.sources[[str.tiers.end]] <- intervention.matrix(paste0("social", last.forecast.date, ".csv"))
    matrix.used[[str.tiers.end]] <- intervention.matrix.out("social_mat")
}
if(!is.null(lockdown.end)){
    use.previous.matrix(str.tiers.begin, str.holidays.start)
    matrix.series(last.forecast.date, str.tiers.begin, holidays.end - 1, last.forecast.week)
} else {
    matrix.series(last.forecast.date, str.tiers.begin, start.date + ndays - 1, last.forecast.week)
}
if(!is.null(holidays.end)){
    ## matrix.sources[[str.holidays.end]] <- intervention.matrix(paste0("current", str.holidays.end, ".csv"))
    ## matrix.used[[str.holidays.end]] <- intervention.matrix.out(paste0("current_mat", length(str.periodic.breaks) + 2))
    use.previous.matrix(str.holidays.end, str.lockdown.end)
}

## Move from list back to model spec
stopifnot(length(matrix.sources) == length(matrix.used) - 1)
### Order correctly
matrix.sources <- matrix.sources[order(ymd(names(matrix.sources)))]
matrix.used <- matrix.used[order(ymd(names(matrix.used)))]
### Parse back into variables
stopifnot(all(names(matrix.used)[-1] == names(matrix.sources)))
cm.breaks <- as.integer(ymd(names(matrix.sources)) - start.date + 1)
cm.lockdown <- as.character(matrix.sources)
cm.bases <- as.character(matrix.used)

## Where to put tmp files and outputs
projections.basedir <- file.path(out.dir, paste0("projections", scenario.num))
projections.file <- paste0("projections_", scenario.name, ".RData")

## Use the below to add extra beta values if required
#extra.beta <- ymd(20201202) - start.date + 1
#extra.beta.values <- log(beta.changes[[regions]])
extra.beta.values <- extra.beta <- NULL
stopifnot(length(extra.beta) == length(extra.beta.values))
## ## ----------------------------------------------------------

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


## ## mod_pars.Rmd specifications that will change - should only be breakpoints and design matrices
## if(prev.flag & all(prior.r1 == 1))
    value.r1 <- 7
bank.holiday.days.new <- NULL
## ## ---------------------------------------------------------------------------------------------


if(!file.exists(projections.basedir))
    dir.create(projections.basedir)

## ## ## CHANGES TO VARIABLES BASED ON mod_inputs-LIKE SPECIFICATIONS
prev.flag <- 0
sero.flag <- 0 ## Are we interested in serological outputs? Switched off for the moment.
start.hosp <- 1
start.gp <- 1
start.prev <- 1
end.hosp <- ifelse(hosp.flag, ndays, 1)
end.gp <- ifelse(gp.flag, ndays, 1)
end.prev <- ifelse(prev.flag, ndays, 1)
dths.flag <- hosp.flag
cases.flag <- gp.flag

## Create matrix files if required
stopifnot(length(cm.bases) - 1 == length(cm.lockdown))
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
mult.order <- rep(1, length(cm.bases))
cm.mults <- c(cm.mults,
              file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_mod", contact.model, "levels", mult.order, ".txt"))
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

## Checks
stopifnot(length(cm.bases) == length(cm.mults))
stopifnot(all(!is.na(cm.mults)))
stopifnot(length(cm.bases) == length(cm.breaks) + 1)

## The mod_inputs.txt file wont change with each projections so can render it now
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
if(beta.rw.flag) {
	if (is.null(extra.beta)) {
		symlink.design("beta.design.txt")
	} else if (length(extra.beta) == 1) {
		beta.design <- cbind(beta.design, rep(0, dim(beta.design)[1]))
		beta.design <- rbind(beta.design, rep(1, dim(beta.design)[2]))
		write_tsv(as.data.frame(beta.design), file.path(projections.basedir, "beta.design.txt"), col_names = FALSE)
	} else {
		stop("Only implemented one extra beta")
	}
}
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
