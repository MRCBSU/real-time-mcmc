source('R/output/combine_region_runs.R')
save.image("model_runs/20200405regions_alone/_OVERALL_/mcmc.RData")

create_quantiles <- function(data) {
	return(quantile(data, c(0.01, 0.05, 0.25, 0.75, 0.95, 0.99, 0.5)))
}


