## Incubation period - best working estimate - mean 5.2 (4.1-7.0)
## Use these as simulation parameters for the latent period

int.effect <- c(0.6435671, 0.52, 0.1753693, 0.00, 0.52, 0.24, 1.0)
names(int.effect) <- c("lo", "med", "hi", "total", "variable", "lshtm", "nothing")

## shape.dL <- 35.1
## rate.dL <- 6.76  ## These values give the desired mean with a variance of 0.768
value.dl <-2
## shape.dL <- 13.3
## rate.dL <- 4.16

## Serial interval - best working estimate - mean 7.5 (5.3-19)
## Use these to deduce the mean infectious period 2*(serial - latent)
## shape.dI <- 4.47
## rate.dI <- 0.972
value.dI <- 0.972
pars.dI <- c(1.43, 0.549)

## Exponential growth rate
value.egr <- 0.14
pars.egr <- c(31.36, 224)

## Ascertainment parameters
value.pgp <- 0.1
pars.pgp <- c(2.12, 15.8)

## Infection to fatality ratio
value.ifr <- 0.007
pars.ifr <- c(21.6, 3070) / 4

## Initial seeding
value.nu <- c(-19, -17.7)
pars.nu <- c(-17.5, 1.25)

## GP Overdispersion
value.eta <- 1.0;
pars.eta <- c(1.0, 0.2);

## Hosp Overdispersion
value.eta.h <- 1.0
pars.eta.h <- c(1.0, 0.2)

## Delay to death
ddelay.mean <- 17.8
ddelay.sd <- 8.9

## Delays in the line listing data
ldelay.mean <- 5.22
ldelay.sd <- 3.59

## Contact model
which.var <- which(sapply(sapply(names(int.effect), function(x) grep(x, scenario.name, fixed = TRUE)), length) > 0)
contact.reduction <- c(1, int.effect[which.var])
contact.dist <- c(1, 3) ## rep(1, 2) if fixed
contact.pars <- c(20.67, 31.00)
