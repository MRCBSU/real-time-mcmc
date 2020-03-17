gp.data <- "../../data/Linelist/linelist20200313.txt"
gp.denom <- "../../data/Linelist/ll_denom20200313.txt"
start.gp <- 15
ndays.gp <- 26

hosp.flag <- 0
hosp.data <- "NULL"
ndays.hosp <- 1

viro.data <- NULL
viro.denom <- NULL

## ## Number of regions (country-level regions in capitals)
regions <- "UNITED KINGDOM"
## ## Vector of age-group descriptions
age.grps <- "All";
## ## Number of days, including lead-in time, analysis of data and short-term projection
ndays <- 39
## Timing of changes to the contact pattern
cm.breaks <- 6
cm.bases <- rep("single_age.txt", 2)
cm.mults <- c(paste0("single_age_mult", 0:1, ".txt"))


############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############

## Get the number of age groups and regions
nages <- length(age.grps)
nregs <- length(regions)

## Make the output directory if necessary
out.dir <- paste0("model_runs/", out.dir)
flg.createfile <- !file.exists(out.dir)
if(flg.createfile) system(paste("mkdir", out.dir))

## Get the population sizes
pop.input <- c(67,000,000)
## source("get_popn.R")
## Remove spaces from region name.
regions <- gsub(" ", "_", regions, fixed = TRUE)

## Contact Model
if(!exists("cm.breaks")) cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
cm.bases <- paste0("../../contact_mats/", cm.bases)
cm.mults <- paste0("../../contact_mats/", cm.mults)
