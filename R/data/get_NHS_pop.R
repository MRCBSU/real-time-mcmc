## Get the population sizes
require(readr)
require(tidyr)
require(dplyr)
dir.data <- "data"
load(build.data.filepath("population", "pop_nhs.RData"))
get.nhs.region <- function(reg, rlist = nhs.regions){
    if(reg %in% names(nhs.regions)){
        return(reg)
    } else if(toupper(reg) %in% names(nhs.regions)) return(toupper(reg))
}
## Check that regions have population specified
for (region in regions) {
    if (!get.nhs.region(region) %in% names(nhs.regions)) {
        stop(paste(region, "is not specified in `nhs.regions`. Options are:",
                   paste0(names(nhs.regions), collapse=", ")))
    }
}
pop.input <- NULL
pdf.all <- NULL
for(reg in regions){
    reg.nhs <- get.nhs.region(reg)
	if (reg == "Scotland" && age.labs[1] == "All") {
		pop.input <- c(pop.input, 5438100)
	} else {
            pop.full <- pop[pop$Name %in% nhs.regions[[get.nhs.region(reg)]] & !is.na(pop$Name), -(1:3), drop = FALSE]
            pop.full <- apply(pop.full, 2, sum)
            if(age.labs[1] == "All"){
                pop.input <- c(pop.input, pop.full["All ages"])
            } else {
                pdf <- data.frame(age = as.numeric(names(pop.full)[-1]), count = pop.full[-1])
                pdf <- pdf %>%
                    mutate(age.grp = cut(pdf$age, age.agg, age.labs, right = FALSE, ordered_result = T)) %>%
                    group_by(age.grp) %>%
                    summarise(count = sum(count))
                pdf.all <- pdf.all %>% bind_rows(pdf %>% mutate(region = reg))
                pop.input <- c(pop.input, pdf$count)
            }
        }
}
