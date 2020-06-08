value.dl <- 2
## shape.dL <- 13.3
## rate.dL <- 4.16

## Serial interval - best working estimate - mean 7.5 (5.3-19)
## Use these to deduce the mean infectious period 2*(serial - latent)
## shape.dI <- 4.47
## rate.dI <- 0.972
value.dI <- 0.972
pars.dI <- c(1.43, 0.549)

## Exponential growth rate
value.egr <- 0.25
pars.egr <- c(31.36, 224)

## Infection to fatality ratio
if(nA > 1){
    value.ifr <- jitter(rep(0.001, nA - 1))
    var.ifr <- rep(0.005, nA - 1)
} else {
    value.ifr <- 0.007
    var.ifr <- 0.005
}
if(nA == 1){
    pars.ifr <-  c(21.6, 3070) / 4  ## c(4.35, 770)
} else {
    means.ifr <- c(1.61e-5, 4.28e-5, 1.89e-4, 9.02e-4, 8.20e-3, 3.11e-2, 6.04e-2) * 11599/15411
    pars.ifr <- as.vector(rbind(rep(1, length(means.ifr)), (1 - means.ifr) / means.ifr))
    pars.ifr[13] <- 9.72
    pars.ifr[14] <- 156
}

## Infection hospitalisation rate
pars.ihr <- c(1, 1)

## Initial seeding
value.nu <- c(-19, -17.7)
pars.nu <- c(-17.5, 1.25)

## GP Overdispersion
value.eta <- 1.0;
pars.eta <- c(1.0, 0.2);

## Hosp Overdispersion
value.eta.h <- 0.3
pars.eta.h <- c(1.0, 0.2)

## Delay to death
ddelay.mean <- 16.0
ddelay.sd <- 12.6

nm <- length(m.param.breakpoints)
contact.dist <- rep(c(1, rep(2, nm)), nr)
contact.pars <- rep(c(4, 4), nr * nm)
contact.reduction <- rep(1, nr * (nm + 1))
contact.proposal <- rep(c(1, rep(0.001, nm)), nr)
