standalone <- T

fns <- "deaths_comparison_prev.RData"
projections <- FALSE

pth <- getwd()
rmarkdown::render(file.path(Rfile.dir, "vacc_lives_saved.Rmd"),
                  output_file = paste0("vacc_lives_saved", ifelse(projections, "_proj", ""), ".html"),
                  output_dir = pth,
                  clean = T)

