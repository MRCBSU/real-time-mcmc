load("tmp.RData")
library(tidyverse)
out.dir <- getwd()
arg.index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(arg.index))
    arg.index <- 1
args <- commandArgs(trailingOnly = TRUE)
projections_file <- args[arg.index]

regions.total.population <- t(matrix(pop.input, nA, length(regions)))
colnames(regions.total.population) <- age.labs
rownames(regions.total.population) <- regions
population <- as_tibble(regions.total.population, rownames = "region") %>%
  pivot_longer(-region, names_to = "age")

require(rmarkdown)

Rfile.loc <- file.path(file.loc, "R/output")
external <- FALSE

output_filename <- gsub("RData", "html", projections_file)
rmarkdown::render(
    file.path(Rfile.loc, 'projinf-report.Rmd'),
    html_document(pandoc_args = "--self-contained"),
    output_dir = out.dir,
    intermediates_dir = file.path(out.dir, paste0(output_filename, "_tmp")),
    output_file = output_filename
)
