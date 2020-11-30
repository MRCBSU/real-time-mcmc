load("tmp.RData")
out.dir <- getwd()
arg.index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(arg.index)) arg.index <- 1
args <- commandArgs(trailingOnly = TRUE)
projections_file <- args[arg.index]

require(rmarkdown)

Rfile.loc <- file.path(file.loc, "R/output")
external <- FALSE

output_filename <- gsub("RData", "html", projections_file)
rmarkdown::render(
	file.path(Rfile.loc, 'projinf-report.Rmd'),
	html_document(pandoc_args = "--self-contained"),
	output_dir = out.dir,
	intermediates_dir = file.path(out.dir, paste0(output_filename, "_tmp"))
)
