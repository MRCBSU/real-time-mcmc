#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################

date.of.runs <- "20200412"    	# What date is in the input file names?
regions <- c("London", "Outside_London")	# Regions under study
data.desc <- "deaths"                          # Needs to be either 'deaths' or 'reports'
# How big an effect should be assumed for the introduction of lockdown?
# Options are "lo", "med", or "high"
scenario.name <- "L_OL_death_variable_tight"

# Choose the name of the subdirectory in model_runs to use
subdir.name <- paste0(date.of.runs, "/", scenario.name)
out.dir <- file.path(proj.dir, "model_runs", subdir.name)	# Value actually used

# Number of days to run the simulation for.
# Including lead-in time, analysis of data and short-term projection
ndays <- 107

#######################################################################
## INPUT SETTINGS
#######################################################################

# The 'gp' stream in the code is linked to confirmed cases data
gp.flag <- 0					# 0 = off, 1 = on
start.gp <- 15					# What day to start running the likelihood on
end.gp <- 15					# Total days of data, or NULL to infer from length of file

# The 'hosp' stream in the code is linked to death data
hosp.flag <- 1					# 0 = off, 1 = on
start.hosp <- 1	 ## 35			# What day to start running the likelihood on
end.hosp <- 51		# Total days of data, or NULL to infer from length of file

viro.data <- NULL
viro.denom <- NULL

# Vector of age-group descriptions
age.grps <- "All"

## CONTACT MATRICES SETTINGS
cm.breaks <- 36							# Day numbers where breaks happen
cm.bases <- rep("single_age.txt", 2)	# Base matrices
# Modifiers (which element of contact_parameters to use)
cm.mults <- c(paste0("single_age_mult", 0:1, ".txt"))



############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############
source(file.path(proj.dir, "R/data/utils.R"))

## Where are the data files?
dir.data <- file.path(proj.dir, "data")

## Map what we call regions (LHS) to the NHS region(s) they contain
## These no longer calculated `on the fly' and should be handled within the data/population folder.
## Use objects nhs.regions and pop
load(build.data.filepath("population", "pop_nhs.RData"))

## Check that regions have ONS population specified
for (region in regions) {
	if (!region %in% names(nhs.regions)) {
		stop(paste(region, "is not specified in `nhs.regions`"))
	}
} 

## Where to store the data outputs.
# If end.gp and/or end.hosp are NULL then read from data files
set.end.date <- function(user.value, data.file) {
	if (is.null(user.value)) {
		return(length(readLines(data.file)))
	} else {
		return(user.value)
	}
}

if(gp.flag){
    gp.data <- build.data.filepath("RTM_format", "linelist", date.of.runs, "_", regions, ".txt")
    gp.denom <- build.data.filepath("RTM_format", "ll_denom", date.of.runs, "_", regions, ".txt")
} else {
    gp.data <- "NULL"
    gp.denom <- "NULL"
}
if(is.null(end.gp))
    end.gp <- ifelse(gp.flag, set.end.date(end.gp, gp.data), start.gp)

if(hosp.flag){
    hosp.data <- build.data.filepath("RTM_format/deaths",
                                     data.desc,
                                     date.of.runs,
                                     "_",
                                     regions,
                                     ".txt")
} else hosp.data <- "NULL"

if(is.null(end.hosp))
    end.hosp <- ifelse(hosp.flag, set.end.date(end.hosp, hosp.data), start.hosp)

## Data file locations: shouldn't need to be changed, calculated based on above
dir.data <- file.path(proj.dir, "data")
## Get the number of age groups and regions
nages <- length(age.grps)
nregs <- length(regions)


## Contact Model
if(!exists("cm.breaks")) cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
cm.bases <- file.path(proj.dir, "contact_mats", cm.bases)
cm.mults <- file.path(proj.dir, "contact_mats", cm.mults)
