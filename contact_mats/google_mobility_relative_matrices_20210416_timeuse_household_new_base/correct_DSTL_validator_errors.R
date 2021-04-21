library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
in_file <- args[1]
out_file <- args[2]

read_csv(in_file, col_types = cols(AgeBand = col_character())) %>%
  replace_na(list(
    AgeBand = "All",
	Scenario = "Nowcast"
  )) %>%
  mutate(
    Model = recode(Model, `Region/age` = "Regional/age"),
    Scenario = recode(Scenario, MRP = "MTP"),
	#ValueType = recode(ValueType,
	  #infections_prev = "prevalence_mtp",
	  #incidence_prev = "prevalence"
	#)
  ) %>%
  filter(ValueType != "infections_cum") %>%
  write_csv(out_file)
