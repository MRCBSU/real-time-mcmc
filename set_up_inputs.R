regions <- c("London", "Outside_London")
gp.data <- "../../data/Linelist/linelist20200319.txt"
gp.denom <- "../../data/Linelist/ll_denom20200319.txt"
hosp.data <- paste0("../../data/deaths/deaths20200319_", regions, ".txt")

gp.flag <- 0
start.gp <- 15
end.gp <- 29

hosp.flag <- 1
start.hosp <- 1
end.hosp <- 28

viro.data <- NULL
viro.denom <- NULL

## ## Number of regions (country-level regions in capitals)
ons.regions <- list("London" = "LONDON",
                    "Outside_London" = c("NORTH EAST", "NORTH WEST", "YORKSHIRE AND THE HUMBER", "EAST MIDLANDS", "WEST MIDLANDS", "EAST", "SOUTH EAST", "SOUTH WEST")
                    )
## ## Vector of age-group descriptions
age.grps <- "All";
## ## Number of days, including lead-in time, analysis of data and short-term projection
ndays <- 85
## Timing of changes to the contact pattern
cm.breaks <- 34
cm.bases <- rep("single_age.txt", 2)
cm.mults <- c(paste0("single_age_mult", 0:1, ".txt"))


############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############

## Get the number of age groups and regions
nages <- length(age.grps)
nregs <- length(regions)

## Make the output directory if necessary
cur.dir <- getwd()
out.dir <- paste0("model_runs/", out.dir)
flg.createfile <- !file.exists(out.dir)
if(flg.createfile) system(paste("mkdir", out.dir))

## Get the population sizes
require(readr)
require(tidyr)
pop <- read_csv("data/population/popn2018_all.csv")
pop.input <- NULL
for(reg in regions){
    pop.full <- pop[pop$Name %in% ons.regions[[reg]] & !is.na(pop$Name), -(1:3), drop = FALSE]
    pop.full <- apply(pop.full, 2, sum)
    if(age.grps == "All")
        pop.input <- c(pop.input, pop.full["All ages"])
    }
## source("get_popn.R")
## Remove spaces from region name.
regions <- gsub(" ", "_", regions, fixed = TRUE)

## Contact Model
if(!exists("cm.breaks")) cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
cm.bases <- paste0("../../contact_mats/", cm.bases)
cm.mults <- paste0("../../contact_mats/", cm.mults)
