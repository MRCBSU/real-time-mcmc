require(tidyverse)
require(lubridate)
require(rmarkdown)
require(knitr)
require(abind)
require(parallel)

thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else if (.Platform$GUI == "RStudio" || Sys.getenv("RSTUDIO") == "1") {
                # We're in RStudio
                return(rstudioapi::getSourceEditorContext()$path)
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}
Rfile.loc <- dirname(thisFile())
source(file.path(Rfile.loc, "drw.R"))
colcode.fl <- "Prev765SeroNHSBT_All_ONScutoff_IFR7bp_11wk2_prev14-0PHE_3dose_new_mprior_matrices_20220607_stable_household_admissions_no_deaths"
oldcode.fl <- "Prev765SeroNHSBT_All_ONScutoff_IFR7bp_11wk2_prev14-0PHE_3dose_new_mprior_matrices_20220607_stable_household_admissions_no_deaths_chain2"

## load("mcmc.RData")
load(file.path(colcode.fl, "tmp.RData"))
out.dir <- target.dir <- getwd()
## results.dir <- out.dir
## source(file.path(Rfile.loc, "tracePlots.R"))
## source("sim_func.R")

QUANTILES <- c(0.025, 0.5, 0.975)
## sanmitra.fl <- "/scratch/sanmitra/covid19/results/Sanmitra_chain.csv"

## ## mod_inputs.Rmd items that will change in the projections.

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
## ## Compile a forecast output
combine.rtm.output <- function(x, strFld){
    oList <- lapply(x, function(x) do.call(abind, args = list(x[[strFld]], along = 3)))
    oList <- do.call(abind, args = list(oList, along = 0))
    
    }

## ## ## --------------------

## if(!file.exists(file.path(out.dir, "projections")))
##     dir.create(file.path(out.dir, "projections"))

## MCMC control
## num.iterations <- 1
## thin.params <- 1
## thin.outputs <- 1
## adaptive.phase <- 0
## burnin <- 0
## num.threads <- 1

## ## Adjust for moved file location
## old.str <- "real-time-mcmc-old"
## new.str <- "real-time-mcmc"
## replace.loc <- function(var)
##     gsub(old.str, new.str, var, fixed = TRUE)
## old.repo <- "/project/pandemic_flu/Wuhan_Coronavirus/real-time-mcmc-old/"
## new.repo <- gsub(old.str, new.str, old.repo)

## cm.bases <- gsub(new.repo, old.repo, cm.bases)
## inputs.template.loc <- replace.loc(inputs.template.loc)
## proj.dir <- replace.loc(proj.dir)
## hosp.data <- replace.loc(hosp.data)
## sero.data <- lapply(sero.data, replace.loc)

## The mod_inputs.txt file wont change with each projections so can render it now
## render(inputs.template.loc, output_dir = file.path(out.dir, "projections"), output_format = plain_document())
## knit(inputs.template.loc, output = file.path(out.dir, "projections", "mod_inputs.txt"))
## ## ## ------------------------------------------------------------

## ## ## CHANGES TO VARIABLES BASED ON mod_pars.txt-LIKE SPECIFICATIONS

## ## Copy to projection directory other design matrices
## if(!exists("single.ifr")) single.ifr <- TRUE
## if(gp.flag)
##     symlink.design("icr.design.txt")
## if(rw.flag)
##     symlink.design("m.design.txt")
## if(beta.rw.flag)
##     symlink.design("beta.design.txt")
## if(!single.ifr)
##     symlink.design("ifr.design.txt")
## if(vacc.flag){
##     symlink.design("vac.alpha1.design.txt")
##     symlink.design("vac.alphan.design.txt")
##     symlink.design("vac.pi1.design.txt")
##     symlink.design("vac.pin.design.txt")
## }

## ## ## --------------------------------------------------------------

## ## ## MAIN PROJECTION LOOP
new.params <- new.env()
load(file.path(colcode.fl, "mcmc.RData"), envir = new.params)
params <- new.params$params

## ## Get number of iterations
niter <- min(sapply(params, nrow))

## ## For each iteration
## paul.lfx <- mclapply(1:niter, sim_rtm, mcmc = params, mc.cores = detectCores() - 1)
## paul.lfx <- unlist(paul.lfx)
paul.lfx <- new.params$lfx

paul.lrw <- mclapply(1:niter, drw, pars = params, nc = 2, mc.cores = detectCores() - 1)
paul.lrw <- unlist(paul.lrw)
## ## ## Load in Sanmitra's posteriors
## sanmitra.chain <- read_csv(sanmitra.fl, col_names = FALSE)
## names(sanmitra.chain) <- c("infectious_period",
##                            "sero_test_sensitivity",
##                            "sero_test_specificity",
##                            "hosp_negbin_overdispersion",
##                            paste0("prop_case_to_hosp", 1:ncol(params$prop_case_to_hosp)),
##                            "log_beta_rw_sd",
##                            "exponential_growth_rate1",
##                            "log_p_lambda_01",
##                            paste0("contact_parameters", 2:4),
##                            paste0("log_beta_rw", 2:17),
##                            "exponential_growth_rate2",
##                            "log_p_lambda_02",
##                            paste0("contact_parameters", 6:8),
##                            paste0("log_beta_rw", 19:34),
##                            "exponential_growth_rate3",
##                            "log_p_lambda_03",
##                            paste0("contact_parameters", 10:12),
##                            paste0("log_beta_rw", 36:51),
##                            "exponential_growth_rate4",
##                            "log_p_lambda_04",
##                            paste0("contact_parameters", 14:16),
##                            paste0("log_beta_rw", 53:68),
##                            "exponential_growth_rate5",
##                            "log_p_lambda_05",
##                            paste0("contact_parameters", 18:20),
##                            paste0("log_beta_rw", 70:85),
##                            "exponential_growth_rate6",
##                            "log_p_lambda_06",
##                            paste0("contact_parameters", 22:24),
##                            paste0("log_beta_rw", 87:102),
##                            "exponential_growth_rate7",
##                            "log_p_lambda_07",
##                            paste0("contact_parameters", 26:28),
##                            paste0("log_beta_rw", 104:119))
## sanmitra.chain <- sanmitra.chain %>%
##     mutate(contact_parameters1 = 0,
##            contact_parameters5 = 0,
##            contact_parameters9 = 0,
##            contact_parameters13 = 0,
##            contact_parameters17 = 0,
##            contact_parameters21 = 0,
##            contact_parameters25 = 0,
##            log_beta_rw1 = 0,
##            log_beta_rw18 = 0,
##            log_beta_rw35 = 0,
##            log_beta_rw52 = 0,
##            log_beta_rw69 = 0,
##            log_beta_rw86 = 0,
##            log_beta_rw103 = 0)

## sg.params <- list(hosp_negbin_overdispersion = sanmitra.chain %>% select(starts_with("hosp_negbin_overdispersion")),
##                   infectious_period = sanmitra.chain %>% select(starts_with("infectious_period")),
##                   contact_parameters = sanmitra.chain %>% select(starts_with("contact_parameters")),
##                   exponential_growth_rate = sanmitra.chain %>% select(starts_with("exponential_growth_rate")),
##                   log_p_lambda_0 = sanmitra.chain %>% select(starts_with("log_p_lambda_0")),
##                   prop_case_to_hosp = sanmitra.chain %>% select(starts_with("prop_case_to_hosp")),
##                   sero_test_sensitivity = sanmitra.chain %>% select(starts_with("sero_test_sensitivity")),
##                   sero_test_specificity = sanmitra.chain %>% select(starts_with("sero_test_specificity")),
##                   log_beta_rw = sanmitra.chain %>% select(-ends_with("sd")) %>% select(starts_with("log_beta_rw")),
##                   log_beta_rw_sd = sanmitra.chain %>% select(starts_with("log_beta_rw_sd"))
##                   )

## tmp.names <- names(sg.params)
## sg.params <- sapply(names(sg.params), function(x){ if(ncol(sg.params[[x]])>1){
##                                                        return(sg.params[[x]] %>% select(paste0(x, 1:ncol(.))))
##                                                        } else { return(sg.params[[x]]) } }, USE.NAMES = TRUE)

old.params <- new.env()
load(file.path(oldcode.fl, "mcmc.RData"), envir = old.params)
oldcode.params <- old.params$params

## oldcode.lfx <- mclapply(1:niter, sim_rtm, mcmc = oldcode.params, mc.cores = detectCores() - 1)
## oldcode.lfx <- unlist(oldcode.lfx)
oldcode.lfx <- old.params$lfx

oldcode.lrw <- mclapply(1:niter, drw, pars = oldcode.params, nc = 2, mc.cores = detectCores() - 1)
oldcode.lrw <- unlist(oldcode.lrw)

## sg.lfx <- mclapply(1:nrow(sanmitra.chain), sim_rtm, mcmc = sg.params, proj.dir = old.repo, mc.cores = detectCores() - 1)
## sg.lfx <- unlist(sg.lfx)

## sg.lrw <- mclapply(1:nrow(sanmitra.chain), drw, pars = sg.params, nc = 1, mc.cores = detectCores() - 1)
## sg.lrw <- unlist(sg.lrw)

## sg.total.lfx <- sg.lfx + sg.lrw
paul.total.lfx <- paul.lfx + paul.lrw
oldcode.total.lfx <- oldcode.lfx + oldcode.lrw

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
dim.list <- list(region = regions,
                 age = age.labs,
                 date = start.date + 0:(ndays - 1),
                 iteration = seq(2, niter, by = 2)
                 )
## infections <- melt.list(NNI);rm(NNI)
## save.list <- "infections"
## dimnames(infections) <- dim.list

## if(hosp.flag) {
##     deaths <- melt.list(Deaths);rm(Deaths)
##     save.list <- c(save.list, "deaths")
##     dimnames(deaths) <- dim.list
##     }
## if(gp.flag){
##     cases <- melt.list(Cases);rm(Cases)
##     save.list <- c(save.list, "cases")
##     dimnames(cases) <- dim.list
## }
## save(list = save.list, file = file.path(sanmitra.dir, "projections.RData"))

require(ggplot2)
## require(hrbrthemes)

data.p <- data.frame(likelihood = as.vector(paul.lfx), rw = paul.lrw, total = as.vector(paul.lfx) + paul.lrw, code = "chain1", iteration = 1:length(paul.lrw))
data.s <- data.frame(likelihood = as.vector(oldcode.lfx), rw = oldcode.lrw, total = as.vector(oldcode.lfx) + oldcode.lrw, code = "chain2", iteration = 1:length(oldcode.lrw))

df.lik.test <- bind_rows(data.p, data.s) %>% pivot_longer(-c(iteration, code), names_to = "type")


glfx <- df.lik.test %>% ## filter(iteration > 3240) %>%
    ggplot(aes(x=value, after_stat(density), fill = code)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = "identity") +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    theme_minimal() +
    labs(fill="") +
    facet_wrap(~type, scales = "free") ##  +
    ## facet_grid(rows=vars(code), cols=vars(type), scales = "free")

ggsave(file.path(colcode.fl, "likelihood_histograms_check.pdf"), glfx, height = 7, width = 12)

## #### DO THEY RECUR BASED ON THE NEW CODE
## ## ## For each iteration
## render(gsub(old.str, new.str, inputs.template.loc),output_dir = file.path(out.dir, "projections"), output_format = plain_document)

## paul.lfx <- list(old = paul.lfx, new = mclapply(1:niter, sim_rtm, mcmc = params, proj.dir = new.repo, mc.cores = detectCores() - 1))

## paul.lfx$new <- unlist(paul.lfx$new)

## paul.lrw <- mclapply(1:niter, drw, pars = params, nc = 2, mc.cores = detectCores() - 1)
## paul.lrw <- unlist(paul.lrw)

## sg.lfx <- list(old = sg.lfx, new = mclapply(1:nrow(sanmitra.chain), sim_rtm, mcmc = sg.params, proj.dir = new.repo, mc.cores = detectCores() - 1))
## sg.lfx$new <- unlist(sg.lfx$new)

## sg.lrw <- list(old = sg.lrw, new = mclapply(1:nrow(sanmitra.chain), drw, pars = sg.params, nc = 1, mc.cores = detectCores() - 1))
## sg.lrw <- unlist(sg.lrw$new)
rwd <- data.frame(density = c(paul.lrw, oldcode.lrw),
                  iter = c(seq(new.params$burnin + new.params$thin.params, new.params$num.iterations, by = new.params$thin.params),
                           seq(old.params$burnin + old.params$thin.params, old.params$num.iterations, by = old.params$thin.params)),
                  run = c(rep("chain1", length(paul.lrw)), rep("chain2", length(oldcode.lrw)))
                  )

gv <- (rwd %>% ## filter(iter > 450000) %>%
    ggplot(aes(x=iter, y=density, color=run, group=run)) +
    geom_line() +
    facet_wrap(~run) +
    theme_minimal() +
    scale_colour_manual(values=c("#69b3a2", "#404080")) +
    ggtitle("Random-walk density comparison")) %>%
    ggsave(filename=file.path(colcode.fl, "densitycheck.pdf"), width = 12, height = 10)

gchain <- df.lik.test %>% filter(type == "likelihood") %>% ## mutate(code = fct_reorder(code, desc(code))) %>% ## for re-ordering the plotting order to inspect overlapping traces.
    ggplot(aes(x=iteration, y=value, group=code, color=code)) +
    geom_line() +
    ## facet_wrap(~code) +
    theme_minimal() +
    scale_colour_manual(values=c("#69b3a2", "#404080")) +
    ggtitle("Likelihood traces")

ggsave(filename=file.path(colcode.fl, "lfx.trace.pdf"), plot=gchain, width = 12, height = 10)
    
