require(tidyverse)
require(lubridate)
require(rmarkdown)
require(abind)

load("tmp.RData")
load("mcmc.RData")

QUANTILES <- c(0.025, 0.5, 0.975)

## ## mod_inputs.Rmd items that will change in the projections.

## Forecast projection
nforecast.weeks <- 8


## ## ----------------------------------------------------------

## ## mod_pars.Rmd specifications that will change - should only be breakpoints and design matrices

bank.holiday.days.new <- NULL
## ## ---------------------------------------------------------------------------------------------

## ## ## FUNCTION DEFINITIONS
repeat.last.row <- function(real.fl, dummy.fl){
    tmpdata <- read_tsv(real.fl, col_names = FALSE)
    dummy.fl <- file.path(out.dir, "projections", paste0(dummy.fl, ".txt"))
    tmpdata <- bind_rows(tmpdata,
                         tmpdata[rep(nrow(tmpdata), ndays - nrow(tmpdata)), ]) %>%
            mutate(Date = start.date - 1 + (1:ndays)) %>%
            select(-X1) %>%
            select(Date, everything()) %>%
            write_tsv(dummy.fl, col_names = FALSE)
    return(dummy.fl)
}
symlink.design <- function(design)
    file.symlink(file.path(out.dir, design), "projections")
## ## ## --------------------

if(!file.exists(file.path(out.dir, "projections")))
    dir.create(file.path(out.dir, "projections"))

## ## ## CHANGES TO VARIABLES BASED ON mod_inputs-LIKE SPECIFICATIONS
sero.flag <- 0 ## Are we interested in serological outputs? Switched off for the moment.
if (!exists("gp.flag")) gp.flag <- 0
if (!exists("dths.flag")) dths.flag <- hosp.flag
if (!exists("cases.flag")) cases.flag <- gp.flag
ndays <- lubridate::as_date(date.data) - start.date + (7 * nforecast.weeks) + 1
start.hosp <- 1
start.gp <- 1
end.hosp <- ifelse(hosp.flag, ndays, 1)
end.gp <- ifelse(gp.flag, ndays, 1)

cm.breaks <- append(cm.breaks, 249)
cm.bases <- append(cm.bases, cm.bases[length(cm.bases)])
cm.mults <- append(cm.mults, "/home/joshuab/real-time-mcmc/contact_mats/ag1_mult2.txt")

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

## MCMC control
num.iterations <- 1
thin.params <- 1
thin.outputs <- 1
adaptive.phase <- 0
burnin <- 0
num.threads <- 1

## The mod_inputs.txt file wont change with each projections so can render it now
render(inputs.template.loc, output_dir = file.path(out.dir, "projections"), output_format = plain_document)

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
        write_tsv(file.path(out.dir, "projections", "d_o_w_design_file.txt"), col_names = FALSE)
}

## Copy to projection directory other design matrices
if(gp.flag)
    symlink.design("icr.design.txt")
if(rw.flag)
    symlink.design("m.design.txt")
if(beta.rw.flag)
    symlink.design("beta.design.txt")
## ## ## --------------------------------------------------------------

## ## ## MAIN PROJECTION LOOP

## ## Set-up output quantities
NNI <- NNI.files <- vector("list", nr)
if(hosp.flag) Deaths <- Deaths.files <- vector("list", nr)
if(gp.flag) Cases <- Cases.files <- vector("list", nr)

## ## Get number of iterations
niter <- min(sapply(params, nrow))

## ## For each iteration
setwd("projections")
pct <- 0
for(iter in 1:niter){
    
    ## Get the parameters into the right variable names.
    if(gp.flag) value.eta <- params$gp_negbin_overdispersion[iter, , drop = FALSE]
    if(hosp.flag) value.eta.h <- params$hosp_negbin_overdispersion[iter, , drop = FALSE]
    value.dI <- params$infectious_period[iter, , drop = FALSE]
    if(gp.flag) value.pgp  <- params$prop_case_to_GP_consultation[iter, , drop = FALSE]
    contact.reduction <- params$contact_parameters[iter, , drop = FALSE]
    beta.rw.vals <- params$log_beta_rw[iter, , drop = FALSE]
	contact.reduction[3] <- log(0.6) - log(posterior.R0[iter, 1]) - sum(beta.rw.vals)
    log_beta_rw_sd <- params$log_beta_rw_sd[iter, -1, drop = FALSE]
    value.egr <- params$exponential_growth_rate[iter, , drop = FALSE]
    value.nu <- params$log_p_lambda_0[iter, , drop = FALSE]
    if(hosp.flag) value.ifr <- params$prop_case_to_hosp[iter, , drop = FALSE]
    if(gp.flag) pars.dow <- params$day_of_week_effects[iter, , drop = FALSE]
    sero.sens <- params$sero_test_sensitivity[iter, , drop = FALSE]
    sero.spec <- params$sero_test_specificity[iter, , drop = FALSE]
    
    ## Compile the mod_pars.txt file
    render(gsub("mod_pars", "sim_mod_pars", pars.template.loc), output_file = "mod_pars.txt",
           output_dir = file.path(out.dir, "projections"), output_format = plain_document, quiet = TRUE)
    
    ## Run the code
    system(file.path(proj.dir, "rtm_optim_1ag"), intern = TRUE)

    ## Read the outputs in and append to output objects
    for(intr in 1:nr)
    {
        NNI.files[[intr]] <- file(paste0("NNI_", regions[intr]), "rb")
        if(dths.flag)
            Deaths.files[[intr]] <- file(paste0("Hosp_", regions[intr]), "rb")
        if(cases.flag)
            Cases.files[[intr]] <- file(paste0("GP_", regions[intr]), "rb")
    }
    names(NNI.files) <- regions
    if(dths.flag) names(Deaths.files) <- regions
    if(cases.flag) names(Cases.files) <- regions

    for(intr in 1:r)
    {
        NNI[[intr]] <- NNI[[intr]] %>%
            abind(readBin(NNI.files[[intr]], double(), n = num.iterations * ndays * nA) %>%
                  array(dim = c(nA, ndays, num.iterations)),
                  along = 3)
        close(NNI.files[[intr]])
        if(dths.flag){
            Deaths[[intr]] <- Deaths[[intr]] %>%
                abind(readBin(Deaths.files[[intr]], double(), n = num.iterations * ndays * nA) %>%
                      array(dim = c(nA, ndays, num.iterations)),
                      along = 3)
            close(Deaths.files[[intr]])
        }
        if(cases.flag){
            Cases[[intr]] <- Cases[[intr]] %>%
                abind(readBin(Cases.files[[intr]], double(), n = num.iterations * ndays * nA) %>%
                      array(dim = c(nA, ndays, num.iterations)),
                      along = 3)
            close(Cases.files[[intr]])
        }
    }
    if((100 * iter / niter) > (pct + 1)){
        pct <- pct + 1
        cat(pct, "% Complete\n")
    }
}
names(NNI) <- regions
if(dths.flag) names(Deaths) <- regions
if(cases.flag) names(Cases) <- regions
## ## ## --------------------

## ## ## SAVE SOME OUTPUTS
dim.list <- list(age = age.labs,
                 date = start.date + 0:(ndays - 1),
                 iteration = 1:niter
                 )
infections <- abind(NNI, along = 0);rm(NNI)
save.list <- "infections"
dimnames(infections)[-1] <- dim.list
names(dimnames(infections)) <- c("region", names(dim.list))
if(hosp.flag) {
    deaths <- abind(Deaths, along = 0);rm(Deaths)
    save.list <- c(save.list, "deaths")
    dimnames(deaths)[-1] <- dim.list
    names(dimnames(deaths)) <- c("region", names(dim.list))
    }
if(gp.flag){
    cases <- abind(Cases, along = 0);rm(Cases)
    save.list <- c(save.list, "cases")
    dimnames(cases)[-1] <- dim.list
    names(dimnames(cases)) <- c("region", names(dim.list))
}
save(list = save.list, file = "../projections.RData")

## ## ## Housekeeping
lapply(hosp.data, file.remove)
