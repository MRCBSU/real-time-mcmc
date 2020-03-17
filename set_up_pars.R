## Incubation period - best working estimate - mean 5.2 (4.1-7.0)
## Use these as simulation parameters for the latent period
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

## Initial seeding
value.nu <- -19
pars.nu <- c(-17.5, 1.25)

## GP Overdispersion
value.eta <- 1.0;
pars.eta <- c(1.0, 0.2);

## Contact model
contact.reduction <- rep(1, 2)
