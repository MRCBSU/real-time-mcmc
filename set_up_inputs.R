require(readr)
require(lubridate)

#######################################################################
## INPUT SETTINGS
#######################################################################

if(gp.flag){
    start.gp <- 15			# What day to start running the likelihood on
    end.gp <- NULL			# Total days of data, or NULL to infer from length of file
} else {
    start.gp <- 1
    end.gp <- 1
}

## The 'hosp' stream in the code is linked to death data
if(!exists("hosp.flag")) hosp.flag <- 1	# 0 = off, 1 = on
reporting.delay <- ifelse(data.desc == "deaths", reporting.delay, 0)
if(hosp.flag){
    start.hosp <- ifelse(data.desc == "reports", 35, 1) ## 35 # Day number on which to start likelihood calculation
    ## Total days of data, or NULL to infer from length of file
    end.hosp <- lubridate::as_date(date.data) - reporting.delay - start.date + 1
} else {
    start.hosp <- 1
    end.hosp <- 1
}
## The 'sero' stream in the code
if(!exists("sero.flag")) sero.flag <- 1
if(sero.flag){ ## Need to remove dependency  on rtm.plot as it may not necessarily be defined.
	if(exists("rtm.plot")) {
		start.sero <- min(rtm.plot$date) - start.date + 1
		end.sero <- max(rtm.plot$date) - start.date + 1
	} else {
		warning('Running sero likelihood for whole period')
		start.sero <- 1
		end.sero <- ndays 
	}
}
## The 'viro' stream in the code
viro.data <- NULL
viro.denom <- NULL

# Vector of age-group descriptions
if(!exists("age.labs"))
    age.labs <- "All"


## CONTACT MATRICES SETTINGS
## Load Edwin's base matrices from contactsr
matrix.dir <- file.path(
	proj.dir, "contact_mats",
	paste0("google_mobility_relative_matrices_", google.data.date)
)
cm.breaks <- c(36, 43, 50, 57, 64, 71, 78, 85) ## Day numbers where breaks happen
mat.dates <- start.date + cm.breaks - 1
lst <- readRDS(file.path(matrix.dir, "base_matrices.rds"))
lst$England$all$m <- lst$England$all$m * 1e7
cm.files <- "england_8ag_contact.txt"
for(i in 1:length(cm.breaks))
    cm.files <- c(cm.files, paste0("england_8ag_contact_ldwk", i, "_", google.data.date, ".txt"))
cm.bases <- file.path(proj.dir, "contact_mats", cm.files) ## Base matrices
cm.lockdown.fl <- paste0("England", mat.dates, "all.csv")
cm.lockdown <- file.path(matrix.dir, cm.lockdown.fl)
idx <- 1
if(!all(file.exists(cm.bases))){
    adf <- as.data.frame(lst$England$all$m)
    write_tsv(adf, cm.bases[idx], col_names = FALSE)
    for(fl in cm.lockdown){
        idx <- idx + 1
        mat <- read_csv(fl) * adf 
        write_tsv(mat, cm.bases[idx], col_names = FALSE)
    }
}
## Modifiers (which element of contact_parameters to use)
cm.mults <- file.path(proj.dir, "contact_mats", 
                      paste0("ag", nA, "_mult", 0:8, ".txt"))
mult.order <- c(0, rep(1, length(cm.breaks)))
## mult.order <- 0:length(cm.breaks)
if(!all(file.exists(cm.mults))){
    mult.mat <- lapply(unique(mult.order), function(x) matrix(x, nA, nA))
    for(i in 1:length(mult.mat)) write_tsv(as.data.frame(mult.mat[[i]]),
                                           cm.mults[i],
                                           col_names = FALSE)
    }
cm.mults <- cm.mults[mult.order+1]

## MCMC settings
num.iterations <- 550000
stopifnot(num.iterations < 1e6) # mod_inputs.txt format does not support integers >= one million
burnin <- 50000
adaptive.phase <- burnin / 2
thin.outputs <- 98 	# After how many iterations to output each set of NNI, deaths etc.
thin.params <- 49  # After how many iterations to output each set of parameters
stopifnot(thin.outputs %% thin.params == 0) # Need parameters on iterations we have outputs



############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############
dir.data <- file.path(proj.dir, "data")
if(sys.nframe() <= 4){ ## Check if below source files might have already been loaded in
    source(file.path(proj.dir, "R/data/utils.R"))
    source(file.path(proj.dir, "config.R"))
}

## Map what we call regions (LHS) to the NHS region(s) they contain
## These no longer calculated `on the fly' and should be handled within the data/population folder.
## Use objects nhs.regions and pop
load(build.data.filepath("population", "pop_nhs.RData"))

get.nhs.region <- function(reg, rlist = nhs.regions){
    if(reg %in% names(nhs.regions)){
        return(reg)
    } else if(toupper(reg) %in% names(nhs.regions)) return(toupper(reg))
}
## Check that regions have population specified
for (region in regions) {
    if (!get.nhs.region(region) %in% names(nhs.regions) && region != "Scotland") {
        stop(paste(region, "is not specified in `nhs.regions`. Options are:",
                   paste0(names(nhs.regions), collapse=", ")))
    }
}

# If end.gp and/or end.hosp are none then read from data files
set.end.date <- function(user.value, data.file) {
	if (is.null(user.value)) {
		return(max(length(readLines(data.file[1]))))
	} else {
		return(user.value)
	}
}
# Where are the data files?
dir.data <- file.path(proj.dir, "data")
source(file.path(proj.dir, "R/data/utils.R"))

gp.data <- "NULL"
gp.denom <- "NULL"
if (gp.flag == 1) {
	gp.data <- build.data.filepath("RTM_format", "linelist", date.data, ".txt")
	gp.denom <- build.data.filepath("RTM_format", "ll_denom", date.data, ".txt")
	if(is.null(end.gp)) end.gp <- set.end.date(end.gp, gp.data)
}
hosp.data <- "NULL"
if (hosp.flag == 1) {
    hosp.data <- data.files
    if(!all(file.exists(hosp.data))) {
		print(hosp.data[which(!file.exists(hosp.data))])
        stop("Above hospitalisation data files does not exist")
	}
    if(is.null(end.hosp)) end.hosp <- set.end.date(end.hosp, hosp.data)
}
sero.data <- list(sample = "NULL",
                  positive = "NULL")
if (sero.flag == 1) {
    sero.data <- list(sample = serosam.files, positive = seropos.files)
    if(!all(sapply(sero.data, function(x) all(file.exists(x)))))
        stop("One of the specified serology data files does not exist")
    if(is.null(end.sero)) end.sero <- set.end.date(end.sero, sero.data)
}
## Contact Model
if(!exists("cm.breaks")) {cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
cm.bases <- file.path(proj.dir, "contact_mats", cm.bases)
cm.mults <- file.path(proj.dir, "contact_mats", cm.mults)
}
