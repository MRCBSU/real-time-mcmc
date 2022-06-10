require(tidyverse)
require(cubelyr)
require(lubridate)

load("tmp.RData")
load("output_matrices.RData")

### Function may not be needed - check if already loaded in ###
pw.cut.columns <- function(breakpoints, max.days, intervals.per.day, start = 0, end = max.days){

  cutpoints <- c(0, breakpoints, max.days)
  
  as.numeric(cut(start + (1:((end - start) * intervals.per.day) / intervals.per.day), breaks = cutpoints))
}
### End function ###

### Iteration labelling ###
int_iter <- 0:(num.iterations - 1)
## parameter.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.params)
parameter.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= burnin]
## outputs.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.outputs)
outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]

## Parameterise the serological observation model
if(sero.end.date > sero.end.1stwv){
    idx.flag <- pw.cut.columns(sero.end.1stwv - start.date + 1, sero.end.date - start.date + 1, intervals.per.day = 1)
} else { idx.flag <- rep(1, sero.end.date - start.date + 1) }
if("sero_test_sensitivity" %in% names(params)){
    sens <- params$sero_test_sensitivity[, idx.flag]
    spec <- params$sero_test_specificity[, idx.flag]
    qsens <- round(apply(params$sero_test_sensitivity, 2, quantile, probs = c(0.025, 0.5, 0.975)), digits = 3)
    qspec <- round(apply(params$sero_test_specificity, 2, quantile, probs = c(0.025, 0.5, 0.975)), digits = 3)
} else {
    sens <- t(array(sero.sens, dim = c(max(idx.flag), length(parameter.iterations))))[, idx.flag]
    spec <- t(array(sero.spec, dim = c(max(idx.flag), length(parameter.iterations))))[, idx.flag]
    qsens <- round(apply(sens, 2, quantile, probs = c(0.025, 0.5, 0.975)), digits = 3)
    qspec <- round(apply(spec, 2, quantile, probs = c(0.025, 0.5, 0.975)), digits = 3)
}
dimnames(sens) <- dimnames(spec) <- list(iteration = parameter.iterations, date = start.date:sero.end.date)
sens_spec <- sens %>%
    as.tbl_cube(met_name = "sens") %>% as_tibble() %>% mutate(date = as_date(date)) %>% inner_join(spec %>%
                                                                                                   as.tbl_cube(met_name = "spec") %>% as_tibble() %>% mutate(date = as_date(date)))

# write_rds(qsens, "qsens_seeroprev.rds", "xz", compression = 9L)
# write_rds(qspec, "qspec_seeroprev.rds", "xz", compression = 9L)
# stop()
## get the data - filtering only to times of non-zero samples
sero.dat <- left_join(rtm.sam, rtm.pos, by = c("date", "region", "age.grp")) %>%
    rename(n = n.x,
           y = n.y,
           age = age.grp) %>%
    filter(n > 0)

## Calculate the seropositivities
df_infections <- cum_infections %>%
    as.tbl_cube(met_name = "infections") %>%
    as_tibble() %>%
    filter(age %in% unique(sero.dat$age)) %>%
    left_join(population) %>%
    mutate(date = lubridate::as_date(date)) %>%
    rename(population = value) %>%
    left_join(sens_spec)

df_infections <- df_infections %>%
    mutate(pos = infections / population, p = (sens * pos) + ((1 - spec) * (1 - pos))) %>%
    left_join(sero.dat)
tempx <- df_infections %>%
    filter(!is.na(n)) %>%
    mutate(samp = rbinom(n = sum(!is.na(n)), size = n, prob = p))
df_infections <- df_infections %>% left_join(tempx)
rm(tempx)
## Then run the aggregations by age and region
df_infections_age <- df_infections %>%
    group_by(date, age, iteration, sens, spec) %>%
    summarise(infections = sum(infections), population = sum(population), n = sum(n, na.rm = TRUE), y = sum(y, na.rm = TRUE), samp = sum(samp, na.rm = TRUE)) %>%
    mutate(pos = infections / population, obs = y / n) %>%
    group_by(date, age, n, y) %>%
    summarise(pos_lo = quantile(pos, 0.025),
              pos_med = median(pos),
              pos_hi = quantile(pos, 0.975),
              samp_lo = quantile(samp, 0.025),
              samp_med = median(samp),
              samp_hi = quantile(samp, 0.975))

## By region
df_infections_region = df_infections %>%
    group_by(date, region, iteration, sens, spec) %>%
    summarise(infections = sum(infections), population = sum(population), n = sum(n, na.rm = TRUE), y = sum(y, na.rm = TRUE), samp = sum(samp, na.rm = TRUE)) %>%
    mutate(pos = infections / population, obs = y / n) %>%
    group_by(date, region, n, y) %>%
    summarise(pos_lo = quantile(pos, 0.025),
              pos_med = median(pos),
              pos_hi = quantile(pos, 0.975),
              samp_lo = quantile(samp, 0.025),
              samp_med = median(samp),
              samp_hi = quantile(samp, 0.975))

## ## Plotting

write_rds(df_infections_region, "seroprev_plot.rds", "xz", compression = 9L)

## By age
df_infections_early_age <- df_infections_age %>% filter(date < sero.end.1stwv)

gg_age <- df_infections_early_age %>% ggplot(aes(x = date, y = pos_med, ymin = pos_lo, ymax = pos_hi, fill = age)) +
    geom_line() +
    geom_ribbon(alpha=0.5) +
    ylab("Seroprevalence") +
    ylim(c(0,1)) +
    facet_wrap(~age) +
    geom_point(aes(x = date, y = prop, size = n / 100, fill = fit), data = df_infections_early_age %>% filter(n > 0) %>% mutate(prop = y/n, fit = (samp_lo <= y) & (y <= samp_hi)), shape = 21) +
    labs(title = "Attack rate and serological goodness-of-fit",
         subtitle = paste("Sensitivity = ", qsens[2, 1], "(", qsens[1, 1], ",", qsens[3, 1], "); Specificity = ", qspec[2, 1], "(", qspec[1, 1], ", ", qspec[3, 1], ")")) +
    geom_segment(aes(x = date, xend = date, y = samp_lo, yend = samp_hi),
                 data = df_infections_early_age %>% filter(n > 0) %>% mutate(samp_lo = samp_lo / n, samp_hi = samp_hi / n),
                 colour = "grey")
ggsave("gof_byage_early.pdf", gg_age, height = 20, width = 25)

df_infections_late_age <- df_infections_age %>% filter(date >= sero.end.1stwv)

gg_age2 <- df_infections_late_age %>% ggplot(aes(x = date, y = pos_med, ymin = pos_lo, ymax = pos_hi, fill = age)) +
    geom_line() +
    geom_ribbon(alpha=0.5) +
    ylab("Seroprevalence") +
    ylim(c(0,1)) +
    facet_wrap(~age) +
    geom_point(aes(x = date, y = prop, size = n / 100, fill = fit), data = df_infections_late_age %>% filter(n > 0) %>% mutate(prop = y/n, fit = (samp_lo <= y) & (y <= samp_hi)), shape = 21) +
    labs(title = "Attack rate and serological goodness-of-fit",
         subtitle = paste("Sensitivity = ", qsens[2, 1], "(", qsens[1, 1], ",", qsens[3, 1], "); Specificity = ", qspec[2, 1], "(", qspec[1, 1], ", ", qspec[3, 1], ")")) +
    geom_segment(aes(x = date, xend = date, y = samp_lo, yend = samp_hi),
                 data = df_infections_late_age %>% filter(n > 0) %>% mutate(samp_lo = samp_lo / n, samp_hi = samp_hi / n),
                 colour = "grey")

ggsave("gof_byage_late.pdf", gg_age2, height = 20, width = 25)

## By region
df_infections_early_region <- df_infections_region %>% filter(date < sero.end.1stwv)

gg_region <- df_infections_early_region %>% ggplot(aes(x = date, y = pos_med, ymin = pos_lo, ymax = pos_hi, fill = region)) +
    geom_line() +
    geom_ribbon(alpha=0.5) +
    ylab("Seroprevalence") +
    ylim(c(0,1)) +
    facet_wrap(~region) +
    geom_point(aes(x = date, y = prop, size = n / 100, fill = fit), data = df_infections_early_region %>% filter(n > 0) %>% mutate(prop = y/n, fit = (samp_lo <= y) & (y <= samp_hi)), shape = 21) +
    labs(title = "Attack rate and serological goodness-of-fit",
         subtitle = paste("Sensitivity = ", qsens[2, 1], "(", qsens[1, 1], ",", qsens[3, 1], "); Specificity = ", qspec[2, 1], "(", qspec[1, 1], ", ", qspec[3, 1], ")")) +
    geom_segment(aes(x = date, xend = date, y = samp_lo, yend = samp_hi),
                 data = df_infections_early_region %>% filter(n > 0) %>% mutate(samp_lo = samp_lo / n, samp_hi = samp_hi / n),
                 colour = "grey")
ggsave("gof_byregion_early.pdf", gg_region, height = 20, width = 25)

df_infections_late_region <- df_infections_region %>% filter(date >= sero.end.1stwv)

gg_region2 <- df_infections_late_region %>% ggplot(aes(x = date, y = pos_med, ymin = pos_lo, ymax = pos_hi, fill = region)) +
    geom_line() +
    geom_ribbon(alpha=0.5) +
    ylab("Seroprevalence") +
    ylim(c(0,1)) +
    facet_wrap(~region) +
    geom_point(aes(x = date, y = prop, size = n / 100, fill = fit), data = df_infections_late_region %>% filter(n > 0) %>% mutate(prop = y/n, fit = (samp_lo <= y) & (y <= samp_hi)), shape = 21) +
    labs(title = "Attack rate and serological goodness-of-fit",
         subtitle = paste("Sensitivity = ", qsens[2, 1], "(", qsens[1, 1], ",", qsens[3, 1], "); Specificity = ", qspec[2, 1], "(", qspec[1, 1], ", ", qspec[3, 1], ")")) +
    geom_segment(aes(x = date, xend = date, y = samp_lo, yend = samp_hi),
                 data = df_infections_late_region %>% filter(n > 0) %>% mutate(samp_lo = samp_lo / n, samp_hi = samp_hi / n),
                 colour = "grey")

ggsave("gof_byregion_late.pdf", gg_region2, height = 20, width = 25)

## by region -- all time
gg_region_all <- df_infections_region %>% ggplot(aes(x = date, y = pos_med, ymin = pos_lo, ymax = pos_hi)) +
    geom_line() +
    geom_ribbon(alpha = 0.5) +
    ylab("Seroprevalence") +
    ylim(c(0, 1)) +
    facet_wrap(~region) +
    geom_point(aes(x = date, y = prop, size = n / 100, fill = fit), data = df_infections_region %>% filter(n > 0) %>% mutate(prop = y/n, fit = as.factor(ifelse(samp_lo <= y, 1, 0) + ifelse(samp_hi <= y, 1, 0))), shape = 21) +
    labs(title = "Attack rate and serological goodness-of-fit",
         subtitle = paste("Sensitivity = ", qsens[2, 1], "(", qsens[1, 1], ",", qsens[3, 1], "); Specificity = ", qspec[2, 1], "(", qspec[1, 1], ", ", qspec[3, 1], ")")) +
    geom_segment(aes(x = date, xend = date, y = samp_lo, yend = samp_hi),
                 data = df_infections_region %>% filter(n > 0) %>% mutate(samp_lo = samp_lo / n, samp_hi = samp_hi / n),
                 colour = "grey")

ggsave("gof_byregion_all.pdf", gg_region_all, height = 20, width = 25)
