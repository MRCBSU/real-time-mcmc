standalone <- T

fns <- "deaths_comparison_prev.RData"

pth <- getwd()
rmarkdown::render(file.path(Rfile.dir, "vacc_lives_saved.Rmd"),
                  output_dir = pth,
                  clean = T)

