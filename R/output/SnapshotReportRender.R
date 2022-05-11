library(here)

# load the temp rdata file
load("tmp.RData")

# Location the script is called from
here <- here()

# Snapshot analysis to use
str.date <- "2022-03-04"

## Location of this script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

file.loc <- dirname(thisFile())

# Flag to determine version of report to generate
map_flag <- F

## Map files
nhs_file <- file.path(file.loc, "..", "..", "maps", "nhs_region_sf.rds")
ons_file <- file.path(file.loc, "..", "..", "maps", "ons_region_sf.rds")

# Find the file which match the date of the snapshot analysis and format of the name
files <- list.files(pattern = paste("^snapshot", str.date, ".RData$", sep = "*"))

# Read in the snapshot data file
load(files)

# Render the markdown document
if(map_flag) {
        rmarkdown::render(file.path(file.loc, "SnapshotReportMaps.Rmd"), output_dir = here, output_file = paste0(str.date, '_snapshot_report.html'), envir = globalenv())
} else {
        rmarkdown::render(file.path(file.loc, "SnapshotReport.Rmd"), output_dir = here, output_file = paste0(str.date, '_snapshot_report.html'), envir = globalenv())
}