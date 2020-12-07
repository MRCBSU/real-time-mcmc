require(tidyverse)
require(lubridate)
require(rmarkdown)
require(abind)
require(parallel)
require(knitr)

load("mcmc.RData")
load("tmp.RData")

source(file.path(proj.dir, "R/output/sim_func.R"))

QUANTILES <- c(0.025, 0.5, 0.975)

## ## mod_inputs.Rmd items that will change in the projections.

## Forecast projection
nforecast.weeks <- 9
job.num <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
stopifnot(job.num >= 1 && job.num <= 6)
christmas <- ifelse(job.num %in% c(1, 3, 5), "low", "high")
if (job.num <= 2) {
	december <- "high"
} else if (job.num <= 4) {
	december <- "medium"
} else {
	deceber <- "low"
}
scenario.name <- paste0(december, "_december_", christmas, "_christmas_")

google.data.date <- ymd(google.data.date)

intervention.matrix <- function(name) {
  ## Use the next line to specify where the new matrices are stored
  intervention.dir <- file.path(dirname(matrix.dir), "scenarios_20201127")
  return(file.path(intervention.dir, name))
}

intervention.matrix.out <- function(name) {
  file.path(proj.dir, "contact_mats", paste0(cm.file.base, "_", name, ".txt"))
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


matrix.sources <- as.list(cm.lockdown)
names(matrix.sources) <- c(cm.breaks) + start.date - 1
matrix.used <- as.list(cm.bases)
names(matrix.used) <- c(1, cm.breaks) + start.date - 1

if (regions == "Northern_Ireland") {
  holiday.matrix.source <- matrix.sources[["2020-12-21"]]
  holiday.matrix.used <- matrix.used[["2020-12-21"]]
  current.matrix.source <- matrix.sources[["2020-12-07"]]
  current.matrix.used <- matrix.used[["2020-12-07"]]
  remove.matrices.after(ymd("2020-11-25"))
  if (december == "high") {
    use.previous.matrix("2020-11-28", "2020-11-02")
    use.previous.matrix("2020-12-11", "2020-10-12")
    use.previous.matrix("2021-01-04", "2020-12-11")
  }
  if (december == "medium") {
    use.previous.matrix("2020-11-28", "2020-11-02")
    matrix.sources[["2020-12-11"]] <- current.matrix.source
    matrix.used[["2020-12-11"]] <- current.matrix.used
    use.previous.matrix("2021-01-04", "2020-12-11")
  }
  if (december == "low") {
    use.previous.matrix("2020-11-28", "2020-11-02")
  }
  matrix.sources[["2020-12-21"]] <- holiday.matrix.source
  matrix.used[["2020-12-21"]] <- holiday.matrix.used
} else if (regions == "Wales") {
  holiday.matrix.source <- matrix.sources[["2020-12-21"]]
  holiday.matrix.used <- matrix.used[["2020-12-21"]]
  remove.matrices.after(ymd("2020-12-03"))
  matrix.sources[["2020-12-21"]] <- holiday.matrix.source
  matrix.used[["2020-12-21"]] <- holiday.matrix.used
  use.previous.matrix("2021-01-04", "2020-11-30")
  beta.dates <- c(ymd(20201204, 20201223, 20201228))
  if (december == "medium") {
    extra.beta <- beta.dates - start.date + 1
    extra.beta.values <- log(c(1.05, 0.95, 1.05))
  }
  if (december == "low") {
    extra.beta <- beta.dates - start.date + 1
    extra.beta.values <- log(c(1.1, 0.9, 1.1))
  }
}
  

# Include Xmas period
if (christmas == "high") {
	matrix.sources[["2020-12-23"]] <- intervention.matrix("social2020-12-21.csv")
	matrix.used["2020-12-23"] <- intervention.matrix.out(paste0("xmas_", christmas))
	use.previous.matrix("2020-12-28", "2020-12-21")
} # else no change


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


# Where to put tmp files and outputs
projections.file <- paste0("projections_", scenario.name, ".RData")
projections.basedir <- file.path(out.dir, paste0("projections", scenario.name))

## Use the below to add extra beta values if required
#extra.beta <- ymd(20201202) - start.date + 1
#extra.beta.values <- log(beta.changes[[regions]])
if (!exists("extra.beta")) extra.beta.values <- extra.beta <- NULL
stopifnot(length(extra.beta) == length(extra.beta.values))
## ## ----------------------------------------------------------


# Create matrix files if required
stopifnot(length(cm.bases) - 1 == length(cm.lockdown))
if(!all(file.exists(cm.bases))){
    idx.miss <- which(!file.exists(cm.bases))
    for(idx in idx.miss){
        mat <- (read_csv(cm.lockdown[idx - 1]) * adf) %>%
            write_tsv(cm.bases[idx], col_names = FALSE)
    }
}
stopifnot(all(file.exists(cm.bases)))

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
dths.flag <- hosp.flag
cases.flag <- gp.flag
## Get the new contract matrix modifiers to use
mult.order <- rep(1, length(cm.bases))
cm.mults <- c(cm.mults,
              file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_mod3levels", mult.order, ".txt"))
              )[1:length(cm.bases)]
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
	} else {
		for (i in 1:length(extra.beta)) {
		  beta.design <- cbind(beta.design, rep(0, dim(beta.design)[1]))
		  beta.design <- rbind(beta.design, rep(1, dim(beta.design)[2]))
		}
		write_tsv(as.data.frame(beta.design), file.path(projections.basedir, "beta.design.txt"), col_names = FALSE)
	}
}
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
    exe <- "optim"
} else exe <- Sys.info()["nodename"]
cat("rtm.exe = ", exe, "\n")
cat("full file path = ", file.path(proj.dir, paste0("rtm_", exe)), "\n")
xtmp <- mclapply(1:niter, sim_rtm, mc.cores = detectCores() - 2, rtm.exe = exe)

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
#lapply(cases.files, file.remove)
#lapply(denoms.files, file.remove)
if(prev.flag)
    lapply(prev.data, file.remove)
rm(list=ls())
