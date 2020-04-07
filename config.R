#######################################################################
## THIS FILE CONTAINS GENERAL PARAMETERS NEEDING TO BE UPDATED
#######################################################################

date.data <- "20200403"


# Number of days to run the simulation for.
# Including lead-in time, analysis of data and short-term projection
ndays <- 62


regions <- "ENGLAND"
# Choose the name of the subdirectory in model_runs to use
subdir.name <- paste0(date.data, "regions_alone")
out.dir <- file.path(proj.dir, "model_runs", subdir.name, regions)	# Value actually used

scenario.name <- "variable"
