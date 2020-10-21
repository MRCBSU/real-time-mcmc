library(tidyverse)

load("tmp.RData")
load("output_matrices.RData")

df_infections <- cum_infections %>%
    as.tbl_cube(met_name = infections) %>%
    as_tibble() %>%
    left_join(population) %>%
    mutate(date = lubridate::as_date(date)) %>%
    rename(population = value)

## get the data - filtering only to times of non-zero samples
sero.dat <- left_join(rtm.sam, rtm.pos, by = c("date", "region", "age.grp")) %>%
    rename(n = n.x,
           y = n.y,
           age = age.grp) %>%
    filter(n > 0)

## only care about the cumulative infections for the relevant region/date combinations



sens_spec <- cbind(params$sero_test_sensitivity, params$sero_test_specificity)[seq(2, nrow(params$sero_test_sensitivity), by = 2), ]
dimnames(sens_spec) <- list(iteration = unique(df_infections$iteration), parameter = c("sens", "spec"))
sens_spec <- sens_spec %>%
    as.data.frame() %>%
    rownames_to_column(var = "iteration") %>%
    mutate(iteration = as.integer(iteration))

df_infections <- df_infections %>%
    left_join(sens_spec) %>%
    mutate(pos = infections / population) %>%
    mutate(p = (sens * pos) + ((1 - spec) * (1 - pos)))
