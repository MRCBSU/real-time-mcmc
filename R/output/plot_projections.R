load("tmp.RData")
out.dir <- getwd()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
	projections_file <- args[1]
}

require(rmarkdown)

Rfile.loc <- file.path(file.loc, "R/output")
external <- FALSE

rmarkdown::render(
	file.path(Rfile.loc, 'projinf-report.Rmd'),
	html_document(pandoc_args = "--self-contained"),
	output_dir = out.dir
)
