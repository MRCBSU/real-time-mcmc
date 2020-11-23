## Get the variables as they were saved on the original computer
setup.env <- new.env()
load("tmp.RData", envir = setup.env)

## Want to change file locations from in.root to out.root
in.root <- "/home/joshuab"
out.root <- "/home/jbb50/rds/hpc-work"

## Get all variable names
var.list <- eapply(setup.env, typeof)
var.names <- names(var.list)[unlist(var.list) == "character"]

## Within the specified environment...
with(setup.env, {
    for(vars in var.names)
        assign(vars, gsub(in.root, out.root, get(vars), fixed = TRUE))
})

## Temporary line to be deleted
## expit <- function(x) exp(x)/(1+exp(x))
## abreaks.icr <- 3:7

save(list = ls(envir = setup.env), file = "tmp.RData", envir = setup.env)
