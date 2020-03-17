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

## Fetch latest C++ files.
comp.root.dir <- "/home/phe.gov.uk/paul.birrell/"
proj.dir <- paste0(comp.root.dir, "Documents/PHE/hdrive/project/pandemic_flu/Real_Time_Modelling/")
c.loc <- paste0(proj.dir, "MCMC/real-time-mcmc/")

## Get the number of age groups and regions
nages <- length(age.grps)
nregs <- length(regions)

## Make the output directory if necessary
out.dir <- paste0("model_runs/", out.dir)
flg.createfile <- !file.exists(out.dir)
if(flg.createfile) system(paste("mkdir", out.dir))
## Populate the directory with the necessary C++ files and compile it
if(flg.createfile){
    system(paste0("ln -s ", c.loc, "*.cc ", out.dir))
    system(paste0("ln -s ", c.loc, "*.h ", out.dir))
    system(paste0("ln -s ", c.loc, "GMakefile ", out.dir))
    ## Change the hard-wiring of the number of age groups.
    header <- readLines(paste0(out.dir, "RTM_Header.h"))
    intHea <- grep("NUM_AGE_GROUPS", header)
    header[intHea] <- paste0("#define NUM_AGE_GROUPS (", nages, ")")
    write(header, file = paste0(out.dir, "RTM_Header.h"), append = F)
}
## And compile the code to get an executable
cur.dir <- getwd()
setwd(out.dir)
system("make -f GMakefile")
setwd(cur.dir)

## Get the population sizes
require(readr)
require(tidyr)
setwd("~/Documents/PHE/stats/Wuhan 2019 Coronavirus/Data/population/")
pop <- read_csv("popn2018_all.csv")
setwd(cur.dir)
pop.input <- NULL
for(reg in regions){
    pop.full <- pop[pop$Name == reg & !is.na(pop$Name), ]
    if(age.grps == "All")
        pop.input <- c(pop.input, pop.full$`All ages`)
    }
## source("get_popn.R")
## Remove spaces from region name.
regions <- gsub(" ", "_", regions, fixed = TRUE)

## Contact Model
if(!exists("cm.breaks")) cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
cm.bases <- paste0("../../contact_mats/", cm.bases)
cm.mults <- paste0("../../contact_mats/", cm.mults)
