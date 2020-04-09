#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################

date.data <- "20200407"


# Number of days to run the simulation for.
# Including lead-in time, analysis of data and short-term projection
ndays <- 75

all.regions <- c(
	"East_of_England",
	"London",
	"Midlands",
	"North_East_and_Yorkshire",
	"North_West",
	"South_East",
	"South_West"
)			# Regions under study
args <- commandArgs(trailingOnly = TRUE)
region.index <- as.integer(args[length(args)])
if (length(region.index) == 0) region.index <- 1
regions <- all.regions[region.index]

# Choose the name of the subdirectory in model_runs to use
subdir.name <- paste0(date.data, "_regions_alone_new")
out.dir <- file.path(proj.dir, "model_runs", subdir.name, regions)	# Value actually used

scenario.name <- "variable"
combined.dir <- file.path(proj.dir, "model_runs", subdir.name, "_OVERALL_")	# Value actually used
