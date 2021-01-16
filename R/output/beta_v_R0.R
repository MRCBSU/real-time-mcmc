library(tidyverse)
library(lubridate)
library(latticeExtra)
library(pammtools)
library(cubelyr)

load("mcmc.RData")

fl.eval <- "dominant_eigenvalues.rds"
if(file.exists(fl.eval))
    Rs <- readRDS("dominant_eigenvalues.rds") %>%
        dplyr::filter(Region == "England",
                      Date >= ymd("20200323"),
                      Date < today(),
                      type == "all") %>%
        mutate(lvalue = log(value))

nw <- nrow(beta.design) / nr
niter <- nrow(params$log_beta_rw)
beta <- exp(beta.design %*% t(params$log_beta_rw)) %>%
    array(dim = c(nw, nr, niter))
names(dim(beta)) <- c("Date", "Region", "iteration")
dimnames(beta) <- list(Date = seq(ymd("20200323"), by = 7, length.out = nw),
                       Region = regions,
                       iteration = 1:niter)

beta <- as.tbl_cube(beta, met_name = "value") %>%
    as_tibble() %>%
    mutate(Date = lubridate::as_date(Date))

qbeta <- beta %>%
    group_by(Date, Region) %>%
    summarise(value = list(quantile(value, probs = c(0.025, 0.5, 0.975)))) %>%
    mutate(quantile = list(paste0("q", names(value[[1]])))) %>%
    unnest(c(quantile, value)) %>%
    pivot_wider(names_from = quantile, values_from = value)

## obj1 <- xyplot(value ~ Date, data = Rs, type = "l", lwd = 2, col="steelblue")
## obj2 <- xyplot(`q50%` ~ Date, data = qbeta, type = "l", lws = 2, col = "#69b3a2")

## doubleYScale(obj1, obj2, add.ylab2 = TRUE, use.style = FALSE)

## Get into a data fram from which we can do a ribbon plot
df_areaStep <- bind_rows(old = qbeta,
                         new = qbeta %>%
                             group_by(Region) %>%
                             mutate(`q50%` = dplyr::lag(`q50%`, default = 1), `q2.5%` = dplyr::lag(`q2.5%`, default = 1), `q97.5%` = dplyr::lag(`q97.5%`, default = 1),
                                   Date = Date - 0.001),
                         current = qbeta %>%
                             filter(Date == max(qbeta$Date)) %>%
                             mutate(Date = ymd(date.data)),
                         .id = "source") %>%
    arrange(Date, source)

require(plotly)
require(htmlwidgets)

gqb <- ggplot(df_areaStep, aes(x = Date, y = `q50%`, colour = Region)) +
    geom_step() +
    geom_ribbon(aes(ymin = `q2.5%`, ymax = `q97.5%`, fill = Region), alpha = 0.4) +
    xlab("Date") +
    ylab("normalised beta")##  +
    ## ylim(c(0, 4))

if(exists("Rs")) gqb <- gqb +
    geom_step(data = Rs %>% mutate(scaled_value = 3 + (value - min(value)) / (max(value) - min(value))), aes(x = Date, y = scaled_value), linetype = 1, lwd = 2, col = "black")
if(!exists("Rs")) gqb <- gqb + facet_wrap(~Region, ncol = 1, scales="free_y") +
                      theme_minimal() +
                      theme(legend.position="none")

ggsave("beta_and_Rstar.pdf", gqb, width = 9, height = 12.5)

pp <- ggplotly(gqb, tooltip = "text")

saveWidget(pp, file="beta_and_Rstar.html")

gscat <- ggplot(qbeta %>% left_join(Rs %>% select(-Region)), aes(x = `q50%`, y = value, label = Date, group = Region)) +
    ## geom_path() +
    geom_point(aes(colour = Date), size = 2.5) +
    ## geom_text()
    theme_classic() +
    ## scale_colour_distiller(palette = "Spectral") +
    xlab("Beta") +
    ylab("R") +
    facet_wrap(~Region, scale = "free_x")

ggsave("beta_Rstar_regional_scatter_nopath.pdf", gscat, width = 12.5, height = 8.25)

gscat_path <- ggplot(qbeta %>% left_join(Rs %>% select(-Region)), aes(x = `q50%`, y = value, label = Date, group = Region)) +
    geom_path() +
    geom_point(aes(colour = Date), size = 2) +
    ## geom_text()
    xlab("Beta") +
    ylab("R") +
    facet_wrap(~Region, scale = "free_x")

ggsave("beta_Rstar_regional_scatter.pdf", gscat_path, width = 12.5, height = 8.25)


#### REPEAT THE ABOVE BUT FOR INCREMENTS
qbeta.inc <- beta %>%
    group_by(Region, iteration) %>%
    mutate(beta_increment = log(value) - log(lag(value, default = 0))) %>%
    filter(beta_increment < Inf) %>%
    group_by(Date, Region) %>%
    summarise(value = list(quantile(beta_increment, probs = c(0.025, 0.5, 0.975)))) %>%
    mutate(quantile = list(paste0("q", names(value[[1]])))) %>%
    unnest(c(quantile, value)) %>%
    pivot_wider(names_from = quantile, values_from = value) %>%
    rename(lincrement = `q50%`,
           lincrement_low = `q2.5%`,
           lincrement_hi = `q97.5%`)

## df_areaStep_inc <- bind_rows(old = qbeta.inc,
##                          new = qbeta.inc %>%
##                              group_by(Region) %>%
##                              mutate(`q50%` = dplyr::lag(`q50%`), `q2.5%` = dplyr::lag(`q2.5%`), `q97.5%` = dplyr::lag(`q97.5%`)),
##                          current = qbeta.inc %>%
##                              filter(Date == max(qbeta.inc$Date)) %>%
##                              mutate(Date = max(Rs$Date)),
##                          .id = "source") %>%
##     arrange(Date, source)

## RsStep_inc <- bind_rows(old = Rs,
##                         new = Rs %>%
##                             group_by(Region) %>%
##                             mutate(`q50$` = dplyr::lag(lvalue

Rs.inc <- Rs %>%
    mutate(lincrement = lvalue - dplyr::lag(lvalue)) %>%
    filter(lincrement < Inf) %>%
    select(-value, -type, -lvalue) %>%
    mutate(lincrement_low = lincrement,
           lincrement_hi = lincrement)

df.incs <- bind_rows(beta = qbeta.inc,
                     R = Rs.inc,
                     .id = "source")

## Lollipop chart
library(hrbrthemes)
extrafont::loadfonts()
glol <- ggplot(data = df.incs %>% mutate(Region = factor(Region)), aes(group = Region)) +
    geom_linerange(aes(x = Date, ymin = 0, ymax=lincrement), colour = "grey", position = position_dodge(width = 4.25)) +
    geom_point(aes(x = Date, y = lincrement, colour = Region), size = 3, position = position_dodge(width = 4.25)) +
    theme_ipsum(base_family = "Helvetica") +
    theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
    ) +
    xlab("Date") +
    ylab("Increment log(val)") +
    facet_wrap(~source, ncol = 1, scale = "free_y")
ggsave("lollipop_increments.pdf", glol, width = 12.5, height = 8.25)




## gqb_inc_b <- ggplot(df_areaStep_inc, aes(x = Date, y = `q50%`, colour = Region)) +
##     geom_step() +
##     geom_ribbon(aes(ymin = `q2.5%`, ymax = `q97.5%`, fill = Region), alpha = 0.4) +
##     xlab("Date") +
##     ylab("log(beta) increments") +

## gqb_inc_a <- ggplot(data = Rs %>%
##                   mutate(lvalue_increment = lvalue - lag(lvalue, default = 0))
##                   mutate(scaled_value = (lvalue - min(lvalue)) / (max(lvalue) - min(lvalue))), aes(x = Date, y = scaled_value), linetype = 1, lwd = 2, col = "black") +
##     geom_step()

gscat_inc <- ggplot(qbeta.inc %>%
                    left_join(Rs.inc %>% select(-Region), by = "Date") %>%
                    rename(beta = lincrement.x,
                           R = lincrement.y),
                    aes(x = beta, y = R, label = Date, group = Region)) +
    geom_point(aes(colour = Date)) +
    xlab("Beta") +
    ylab("R") +
    facet_wrap(~Region, scale = "free_x")

ggsave("beta_Rstar_regional_inc_scatter_nopath.pdf", gscat_inc, width = 12.5, height = 8.25)
