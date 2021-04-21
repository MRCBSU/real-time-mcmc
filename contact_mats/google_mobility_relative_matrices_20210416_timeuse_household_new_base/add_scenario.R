library(tidyverse)
f <- "~/real-time-mcmc/spi-forecasts/date_20201107/PHE_CAM_deaths_ons20201110_medterms.csv"
read_csv(f) %>%
	mutate(Scenario = "MTP") %>%
	write_csv(f)
