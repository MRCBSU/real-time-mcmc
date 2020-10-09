require(tidyverse)
require(ggplot2)
require(lubridate)

## Location of this script
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
############################
R.dir <- dirname(thisFile())
proj.dir <- dirname(dirname(R.dir))
source(file.path(R.dir, "seir_reporting_functions.R"))
out.dir <- file.path(proj.dir, "model_runs", "20200930", "base_varSens_ifr_60cutoff6day_matrices_20200925_deaths")

load(file.path(out.dir, "tmp.RData"))

mcmc.env <- new.env()
load(file.path(out.dir, "mcmc.RData"), envir = mcmc.env)
params <- mcmc.env$params

## ## Values that should be fixed over all iterations.
delta.t <- 0.5
alp <- value.dl + (1 / delta.t)
base.mats <- mult.mats <- list()
for(i in 1:length(cm.bases)){
    base.mats[[i]] <- as.matrix(read_tsv(cm.bases[i], col_names = FALSE))
    mult.mats[[i]] <- as.matrix(read_tsv(cm.mults[i], col_names = FALSE) + 1)
}
mixing.model <- list(base = base.mats, mults = mult.mats, breaks = cm.breaks, max.ind = max(sapply(mult.mats, max)))
M.intervals <- list()
for(i in 1:(length(cm.breaks) + 1)){
    start.range <- ifelse(i == 1, 1, (cm.breaks[i - 1] / delta.t) + 1)
    end.range <- ifelse(i > length(cm.breaks), ndays, cm.breaks[i]) / delta.t
    M.intervals[[i]] <- start.range:end.range
}
mixing.model$intervals <- M.intervals
rm(base.mats, mult.mats, M.intervals)
## infection to confirmed case ratio
if("prop_case_to_GP_consultation" %in% names(params))
    icr.mat <- model.matrix(ex6)

get.prevalence <- function(iter, r){
    
    ## ## Contents of the loop
    x <- lapply(params, function(x) x[iter, ])
    ## nA - num ages, ndays - num days; nr - num regions
    egr <- x$exponential_growth_rate[r]
    aip <- x$infectious_period + (1 / delta.t)
    nu <- x$log_p_lambda_0[r]
    ntimes <- ndays / delta.t
    times <- 1:ntimes
    if("prop_case_to_GP_consultation" %in% names(x)){
        mat.rows <- nrow(icr.mat) / nr
        mat.rows <- ((r - 1)*mat.rows) + 1:mat.rows
        value.pgp <- icr.mat[mat.rows, ] %*% x$prop_case_to_GP_consultation
        }
        
    ## beta model
    beta.block <- length(beta.breaks) + 1
    beta.rows <- ((r - 1) * beta.block) + (1:beta.block)
    beta <- exp(beta.design[beta.rows, ] %*% x$log_beta_rw)
    beta <- beta[pw.cut.columns(beta.breaks, ndays, 1 / delta.t)]
    R.init <- mcmc.env$R0.func(egr, aip, alp)
    Rt <- rep(R.init, ntimes) * beta
    rho <- 2 * delta.t / aip
    sigma <- delta.t * 2 / alp
    alpha <- exp(egr * delta.t) - 1
    I0.tot <- mcmc.env$I0.func(aip, nu, R.init, value.pgp[1], mcmc.env$all.pop[1, r])
    m <- exp(x$contact_parameters[((r-1)*mixing.model$max.ind) + (1:mixing.model$max.ind)])
    mixing.model$M <- lapply(1:length(mixing.model$base),
                             function(mi) mixing.model$base[[mi]] * m[mixing.model$mults[[mi]]])
    popn <- mcmc.env$regions.total.population[r, ]
    mixing.model$idist <- scaled.init.age.distribution(
        ngm.matrix(mixing.model$M[[1]], popn))
    mixing.model$scaledM <- scaled.mixing.matrix(mixing.model, popn)
    I0 <- I0.tot * mixing.model$idist
    
    I1 <- I0 / (1 + (rho / (alpha + rho)))
    I2 <- I1 * (rho / (alpha + rho))
    E2 <- I1 * ((alpha + rho) / sigma)
    E1 <- E2 * ((alpha + sigma) / sigma)
    
    NNI <- p.lambda <- S <- I.1 <- I.2 <- E.1 <- E.2 <- R <- matrix(0, ntimes, nA)
    
    I.1[1, ] <- I1
    I.2[1, ] <- I2
    E.1[1, ] <- E1
    E.2[1, ] <- E2
    
    S[1, ] <- popn
    R[1, ] <- 0
    
    p.beta <- array(0, c(ntimes, nA, nA))
    
    m.max <- min((1:length(mixing.model$intervals))[sapply(1:length(mixing.model$intervals), function(y) ntimes %in% mixing.model$intervals[[y]])])
    
    for(intm in 1:m.max){
        t.max <- min(as.numeric(ntimes), mixing.model$intervals[[intm]][length(mixing.model$intervals[[intm]])])
        t.min <- mixing.model$intervals[[intm]][1]
        for(t in t.min:t.max)
            p.beta[t, , ] <- Rt[t] * mixing.model$scaledM[[intm]] / aip
    }
    
    p.lambda[1, ] <- calc.p.lambda.row(p.beta[1, , ], I.1[1, ], I.2[1, ], delta.t)
    NNI[1, ] <- p.lambda[1, ] * S[1, ]
    
    for(t in 2:nrow(S)){
        S[t, ] <- S[t - 1, ] * (1 - p.lambda[t -1, ])
        E.1[t, ] <- ((1 - ((2 * delta.t) / alp)) * E.1[t - 1, ]) + (p.lambda[t - 1, ] * S[t - 1, ])
        E.2[t, ] <- ((1 - ((2 * delta.t) / alp)) * E.2[t - 1, ]) + ((2 * delta.t / alp) * E.1[t - 1, ])
        I.1[t, ] <- (I.1[t - 1, ] * (1 - (2 * delta.t / aip))) + ((2 * delta.t / alp) * E.2[t - 1, ])
        I.2[t, ] <- (I.2[t - 1, ] * (1 - (2 * delta.t / aip))) + (I.1[t - 1, ] * 2 * delta.t / aip)
        R[t, ] <- R[t - 1, ] + (I.2[t - 1, ] * 2 * delta.t / aip)
        NNI[t, ] <- p.lambda[t - 1, ] * S[t - 1, ]
        p.lambda[t, ] <- calc.p.lambda.row(p.beta[t, , ], I.1[t, ], I.2[t, ], delta.t)
    }

    prev <- (E.1 + E.2 + I.1 + I.2)[seq(from = 2, by = 2, length.out = ndays), ]
    dimnames(prev) <- list(date = ymd("20200216") + days(1:ndays),
                           age = age.labs)

    prev <- prev %>%
        as_tibble() %>%
        mutate(iteration = iter,
               region = regions[r]) %>%
        rownames_to_column("date") %>%
        mutate(date = ymd("20200216") + days(as.numeric(date)))
    
    return(prev)

}

## ## Quantities that are to be looped over, iteration and region.
max.iter <- nrow(params$log_beta_rw)
iters <- sample.int(max.iter, 2000)
prev.array <- expand.grid(iters, 1:nr)
colnames(prev.array) <- c("iteration", "region")

require(parallel)

prev.list <- mclapply(1:nrow(prev.array), function(x) get.prevalence(prev.array[x, 1], prev.array[x, 2]),
                      mc.cores = 4)

prevalence <- do.call(bind_rows, prev.list) %>%
    pivot_longer(cols = -c(date, iteration, region), names_to = "age", values_to = "prevalence")
rm(prev.list)

prevalence <- prevalence %>%
    group_by(date, iteration, age) %>%
    summarise(prevalence = sum(prevalence)) %>%
    mutate(region = "England") %>%
    bind_rows(prevalence)

prev.qts <- prevalence %>%
    group_by(date, age, region) %>%
    summarise(p50 = median(prevalence),
              p2.5 = quantile(prevalence, probs = 0.025),
              p97.5 = quantile(prevalence, probs = 0.975))

save(prev.qts, file = file.path(out.dir, "prevalence.RData"))
## get.prevalence(prev.array[1, 1], prev.array[1, 2])

prev.tot.qts <- prevalence %>%
    group_by(date, region, iteration) %>%
    summarise(prevalence = sum(prevalence)) %>%
    group_by(date, region) %>%
    summarise(p50 = median(prevalence),
              p2.5 = quantile(prevalence, probs = 0.025),
              p97.5 = quantile(prevalence, probs = 0.975)
              )

prev.tot.qts %>%
filter(region == "England",
       date > ymd("20200430")) %>%
ggplot(aes(x = date, y = p50, ymin=p2.5, ymax=p97.5)) +
    geom_line() +
    geom_ribbon(alpha = 0.3)
