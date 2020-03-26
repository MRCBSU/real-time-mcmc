#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################

date.of.runs <- "20200325"    	# What date is in the input file names?
regions <- c("ENGLAND")			# Regions under study

# Choose the name of the subdirectory in model_runs to use
scenario.name <- "med" # Optional component of output directory
subdir.name <- paste0("initial_run_deaths_delaysensENG", date.of.runs, "_", scenario.name)
out.dir <- file.path(proj.dir, "model_runs", subdir.name)	# Value actually used

# Number of days to run the simulation for.
# Including lead-in time, analysis of data and short-term projection
ndays <- 91

#######################################################################
## INPUT SETTINGS
#######################################################################

# The 'gp' stream in the code is linked to confirmed cases data
gp.flag <- 0					# 0 = off, 1 = on
start.gp <- 15					# What day to start running the likelihood on
end.gp <- NULL					# Total days of data, or NULL to infer from length of file

# The 'hosp' stream in the code is linked to death data
hosp.flag <- 1					# 0 = off, 1 = on
start.hosp <- 1					# What day to start running the likelihood on
end.hosp <- NULL				# Total days of data, or NULL to infer from length of file

viro.data <- NULL
viro.denom <- NULL

# Map what we call regions (LHS) to the ONS region(s) they contain
# ONS puts country-level regions in caps
# See the data/popn2018_all.csv file for all possible options
ons.regions <- list(
	"London" = "LONDON",
    "Outside_London" = c("NORTH EAST", "NORTH WEST", "YORKSHIRE AND THE HUMBER",
						 "EAST MIDLANDS", "WEST MIDLANDS", "EAST", "SOUTH EAST",
						 "SOUTH WEST"),
	"UNITED_KINGDOM" = "UNITED KINGDOM",
	"ENGLAND" = "ENGLAND"
)

# Vector of age-group descriptions
age.grps <- "All";

## CONTACT MATRICES SETTINGS
cm.breaks <- 36							# Day numbers where breaks happen
cm.bases <- rep("single_age.txt", 2)	# Base matrices
# Modifiers (which element of contact_parameters to use)
cm.mults <- c(paste0("single_age_mult", 0:1, ".txt"))



############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############

# Check that regions have ONS population specified
for (region in regions) {
	if (!region %in% names(ons.regions)) {
		stop(paste(region, "is not specified in `ons.regions`"))
	}
} 


# Where are the data files?
dir.data <- file.path(proj.dir, "data")
source(file.path(proj.dir, "R/data/utils.R"))
gp.data <- build.data.filepath("RTM_format", "linelist", date.of.runs, ".txt")
gp.denom <- build.data.filepath("RTM_format", "ll_denom", date.of.runs, ".txt")
hosp.data <- build.data.filepath("RTM_format", "deaths", date.of.runs, "_ENGLAND.txt")

# If end.gp and/or end.hosp are none then read from data files
set.end.date <- function(user.value, data.file) {
	if (is.null(user.value)) {
		return(length(readLines(data.file)))
	} else {
		return(user.value)
	}
}
end.gp <- set.end.date(end.gp, gp.data)
end.hosp <- set.end.date(end.hosp, hosp.data)
# Data file locations: shouldn't need to be changed, calculated based on above
dir.data <- file.path(proj.dir, "data")
## Get the number of age groups and regions
nages <- length(age.grps)
nregs <- length(regions)


## Contact Model
if(!exists("cm.breaks")) cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
cm.bases <- file.path(proj.dir, "contact_mats", cm.bases)
cm.mults <- file.path(proj.dir, "contact_mats", cm.mults)
