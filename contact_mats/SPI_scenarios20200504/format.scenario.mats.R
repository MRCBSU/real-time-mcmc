require(readr)
dates <- rep(paste("2020", c("05", "06", "07", "08", "09"), c("11", "01", "01", "15", "04"), sep = "-"), 2)
int.scenario <- c(rep(1, length(dates)/2), rep(2, length(dates)/2))
mat.file <- paste0("scenario", int.scenario, "England", dates, "all.csv")
out.mat.file <- paste0("england_8ag_contact_scenario", int.scenario, "phase", c("1", "2", "3", "4a", "4b"), "_20200501.txt")
lst <- readRDS("base_matrices.rds")
lst$England$all$m <- lst$England$all$m * 1e7

adf <- as.data.frame(lst$England$all$m)
idx <- 0
for(fl in mat.file){
    idx <- which(mat.file %in% fl)
    mat <- read_csv(fl) * adf
    write_tsv(mat, out.mat.file[idx], col_names = FALSE)
}
