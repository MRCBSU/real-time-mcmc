#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################

date.data <- "20200403"

# Choose the name of the subdirectory in model_runs to use
subdir.name <- paste0(date.data, "_regions_deaths")
out.dir <- file.path(proj.dir, "model_runs", subdir.name)	# Value actually used

# Number of days to run the simulation for.
# Including lead-in time, analysis of data and short-term projection
ndays <- 62

regions <- c(
	"East_of_England",
	"London",
	"Midlands",
	"North_East_and_Yorkshire",
	"North_West",
	"South_East",
	"South_West"
)			# Regions under study

