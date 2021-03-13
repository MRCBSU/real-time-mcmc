suppressMessages(extract <- R.utils::extract)
require(tidyverse)

load("mcmc.RData")
load("output_matrices.RData")
load("projections_midterm.RData")

## get the right number of iterations
int_iter <- 0:(num.iterations - 1)
## parameter.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.params)
parameter.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= burnin]
## outputs.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.outputs)
outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]
parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)

## dNNI <- extract(vacc.infections, indices = list(parameter.to.outputs), dims = "iteration") %>%
##     extract(indices = list(1:ndays), dims = "date")
## dNNI <- lapply(regions, function(reg) extract(dNNI, "region" = reg, drop = TRUE) %>%
##                                       aperm(c("age", "date", "iteration")))
## DNNI.files <- vector("list", r)
## for(intr in 1:r){
##     DNNI.files[[intr]] <- file(paste0("Delta_Dis_", regions[intr]), "wb")
##     writeBin(as.vector(dNNI[[intr]]), DNNI.files[[intr]], double())
##     close(DNNI.files[[intr]])
## }

## NNI <- extract(infections, indices = list(parameter.to.outputs), dims = "iteration") %>%
##     extract(indices = list(1:ndays), dims = "date")
## NNI <- lapply(regions, function(reg) extract(NNI, "region" = reg, drop = TRUE) %>%
##                                       aperm(c("age", "date", "iteration")))
## NNI.files <- vector("list", r)
## for(intr in 1:r){
##     NNI.files[[intr]] <- file(paste0("NNI_", regions[intr]), "wb")
##     writeBin(as.vector(NNI[[intr]]), NNI.files[[intr]], double())
##     close(NNI.files[[intr]])
## }

## Deaths <- extract(deaths, indices = list(parameter.to.outputs), dims = "iteration") %>%
##     extract(indices = list(1:ndays), dims = "date")
## Deaths <- lapply(regions, function(reg) extract(Deaths, "region" = reg, drop = TRUE) %>%
##                                       aperm(c("age", "date", "iteration")))
## Deaths.files <- vector("list", r)
## for(intr in 1:r){
##     Deaths.files[[intr]] <- file(paste0("Hosp_", regions[intr]), "wb")
##     writeBin(as.vector(Deaths[[intr]]), Deaths.files[[intr]], double())
##     close(Deaths.files[[intr]])
## }

## Prevalence <- extract(prevalence, indices = list(parameter.to.outputs), dims = "iteration") %>%
##     extract(indices = list(1:ndays), dims = "date")
## Prevalence <- lapply(regions, function(reg) extract(Prevalence, "region" = reg, drop = TRUE) %>%
##                                       aperm(c("age", "date", "iteration")))
## Prevalence.files <- vector("list", r)
## for(intr in 1:r){
##     Prevalence.files[[intr]] <- file(paste0("Prev_", regions[intr]), "wb")
##     writeBin(as.vector(Prevalence[[intr]]), Prevalence.files[[intr]], double())
##     close(Prevalence.files[[intr]])
## }

Sero <- extract(seropos, indices = list(parameter.to.outputs), dims = "iteration") %>%
    extract(indices = list(1:ndays), dims = "date")
Sero <- lapply(regions, function(reg) extract(Sero, "region" = reg, drop = TRUE) %>%
                                      aperm(c("age", "date", "iteration")))
Sero.files <- vector("list", r)
for(intr in 1:r){
    Sero.files[[intr]] <- file(paste0("Sero2_", regions[intr]), "wb")
    writeBin(as.vector(Sero[[intr]]), Sero.files[[intr]], double())
    close(Sero.files[[intr]])
}
