## Get the population sizes
require(readr)
require(tidyr)
require(dplyr)
dir.data <- "data"
load(build.data.filepath("population", "pop_ons.RData"))

regions.lkup <- data.frame(region = as.character(regions), reg.short = toupper(regions), stringsAsFactors = FALSE)
idx <- which(regions.lkup == "East_of_England")
regions.lkup[idx, "reg.short"] <- "EAST"

pop.input <- NULL
for(reg in regions){
	if (reg == "East_of_England") reg <- "EAST"
	pop.full <- pop[pop$Region == toupper(reg),]
	stopifnot(all(pop.full$Age == age.labs))
	stopifnot(nrow(pop.full) == length(age.labs))
	pop.input <- c(pop.input, pop.full$x)
}

pdf.all <- pop %>% left_join(regions.lkup, by = c("Region" = "reg.short")) %>%
    select(-Region) %>%
    rename(age.grp = Age, count = x)
