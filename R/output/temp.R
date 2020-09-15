library(tidyverse)
library(ggplot2)

load("tmp.RData")
load("mcmc.RData")

## Plot the proportion of cases presenting in pillar 2.
age.labs2 <- c("<15yr", age.labs[-(1:3)])
mex6 <- model.matrix(ex6)
pgp <- mex6 %*% t(params$prop_case_to_GP_consultation)
expit <- function(x) exp(x) / (1 + exp(x))
pgp <- expit(pgp)
qpgp <- apply(pgp, 1, quantile, probs = c(0.025, 0.5, 0.975))
qpgp <- array(qpgp, dim = c(dim(qpgp)[1], length(age.labs2), length(regions)))
dimnames(qpgp) <- list(quantile = c("2.5%", "50%", "97.5%"), age = age.labs2, region = regions)
xtmp <- as.tbl_cube(qpgp, met_name = "value") %>%
    as_tibble() %>%
    pivot_wider(id_cols = c(2,3),names_from = quantile, values_from = value)

ggtmp <- ggplot(xtmp) +
    geom_bar(aes(x=region, y=`50%`, group =age),
             colour = "black", stat = "identity", fill = "skyblue", alpha = 0.7, position = "dodge") +
    geom_errorbar(aes(x=region, ymin=`2.5%`, ymax=`97.5%`, group = age),
                  width= 0.5, position = position_dodge(0.9), alpha = 0.9, colour = "orange") +
    ylab("P(infections -> symptomatic pillar 2 case)")

ggsave("case_ascertainment.pdf", ggtmp, height = 9, width = 11)
