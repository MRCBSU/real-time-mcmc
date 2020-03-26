dir.data <- file.path(proj.dir, "data")
source(file.path(proj.dir, "R/data/utils.R"))
gp.data <- build.data.filepath("RTM_format", "linelist", date.of.runs, ".txt")
gp.denom <- build.data.filepath("RTM_format", "ll_denom", date.of.runs, ".txt")
hosp.data <- build.data.filepath("RTM_format", "deaths", date.of.runs, "_ENGLAND.txt")

gp.flag <- 0
start.gp <- 15
end.gp <- 36

hosp.flag <- 1
start.hosp <- 1
end.hosp <- 36

viro.data <- NULL
viro.denom <- NULL

## ## Number of regions (country-level regions in capitals)
ons.regions <- list("London" = "LONDON",
                    "Outside_London" = c("NORTH EAST", "NORTH WEST", "YORKSHIRE AND THE HUMBER", "EAST MIDLANDS", "WEST MIDLANDS", "EAST", "SOUTH EAST", "SOUTH WEST"),
                    "UNITED_KINGDOM" = "UNITED KINGDOM",
                    "ENGLAND" = "ENGLAND"
                    )

## ## Vector of age-group descriptions
age.grps <- "All";
## ## Number of days, including lead-in time, analysis of data and short-term projection
ndays <- 91
## Timing of changes to the contact pattern
cm.breaks <- 36
cm.bases <- rep("single_age.txt", 2)
cm.mults <- c(paste0("single_age_mult", 0:1, ".txt"))


############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############

## Fetch latest C++ files.
c.loc <- file.path(proj.dir, "src")

## Get the number of age groups and regions
nages <- length(age.grps)
nregs <- length(regions)

## Make the output directory if necessary
cur.dir <- getwd()
flg.createfile <- !file.exists(out.dir)
if(flg.createfile) system(paste("mkdir", out.dir))
## Populate the directory with the necessary C++ files and compile it
setwd(out.dir)
if(flg.createfile){
    system(paste0("ln -s ", c.loc, "*.cc ./"))
    system(paste0("ln -s ", c.loc, "*.h ./"))
    system(paste0("ln -s ", c.loc, "GMakefile ./"))
    ## Change the hard-wiring of the number of age groups.
    header <- readLines("RTM_Header.h")
    intHea <- grep("NUM_AGE_GROUPS", header)
    header[intHea] <- paste0("#define NUM_AGE_GROUPS (", nages, ")")
    write(header, file = "RTM_Header.h", append = F)
}
## And compile the code to get an executable
system("make -f GMakefile")
setwd(cur.dir)

## Get the population sizes
require(readr)
require(tidyr)
pop <- read_csv(build.data.filepath("", "popn2018_all.csv"))
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
cm.bases <- file.path(proj.dir, "contact_mats", cm.bases)
cm.mults <- file.path(proj.dir, "contact_mats", cm.mults)
