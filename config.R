#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) region.index <- as.integer(args[length(args)])
if (length(args) > 1) date.data <- args[length(args)-1]


if (!exists("date.data")) date.data <- "20200409"

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
	"South_West",
	"Scotland"
)			# Regions under study
if (!exists("region.index")) region.index <- 1

regions <- all.regions[region.index]
#regions <- "ENGLAND"

# Choose the name of the subdirectory in model_runs to use
data.desc <- "incidence"
scenario.name <- "tight_variable"
subdir.name <- paste0(date.data, "_regions_alone_", scenario.name)
out.dir <- file.path(proj.dir, "model_runs", subdir.name, regions)	# Value actually used
combined.dir <- file.path(proj.dir, "model_runs", subdir.name, "_OVERALL_")	# Value actually used
