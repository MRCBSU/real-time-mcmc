args_in = commandArgs(trailingOnly=TRUE)

print(file.path(args_in[1], "tmp.RData"))

library(tidyverse)

load(file.path(args_in[1], "tmp.RData"))

num.iterations <- 1.7e6L
thin.params <- 200L
thin.outputs <- 400L
adapt.every <- 100L
burnin <- 7e5L

knitr::knit(input = file.path("~/mcmc/real-time-mcmc/inputs/mod_inputs.Rmd"), output = file.path(args_in[1], "mod_inputs_2.txt"))

save.image(file.path(args_in[1], "tmp_2.RData"))