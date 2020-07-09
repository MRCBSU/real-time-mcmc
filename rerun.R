## Run from within output directory
load("tmp.RData")

out.dir <- getwd()

require(rmarkdown)

Rfile.loc <- file.path(file.loc, "R/output")

if (!file.exists("mcmc.RData")) {
	source(file.path(Rfile.loc, "tracePlots.R"))
}
rmarkdown::render(
	file.path(Rfile.loc, 'report-updated.Rmd'),
	html_document(pandoc_args = "--self-contained"),
	output_dir = out.dir,
	clean = FALSE, intermediates_dir = out.dir
)

## Return back to initial directory
setwd(startwd)
