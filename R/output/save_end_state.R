out.dir <- commandArgs(trailingOnly = TRUE)[1]
setwd(out.dir)

load("mcmc.RData")
rm(list = ls()[-grep("params", ls(), fixed = TRUE)])
params <- lapply(params, function(x) x[nrow(x), , drop = FALSE])
save.image("endstate.RData")
