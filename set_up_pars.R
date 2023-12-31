## Incubation period - best working estimate - mean 5.2 (4.1-7.0)
## Use these as simulation parameters for the latent period

int.effect <- c(0.6435671, 0.52, 0.1753693, 0.00, 0.521, 0.24, 1.0)
names(int.effect) <- c("lo", "med", "hi", "total", "variable", "lshtm", "nothing")

## shape.dL <- 35.1
## rate.dL <- 6.76  ## These values give the desired mean with a variance of 0.768
value.dl <-1
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

## Ascertainment parameters
value.pgp <- 0.1
pars.pgp <- c(2.12, 15.8)

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
    means.ifr <- c(1.61e-5, 4.28e-5, 1.89e-4, 9.02e-4, 8.20e-3, 3.11e-2, 6.04e-2)
    pars.ifr <- as.vector(rbind(rep(1, length(means.ifr)), (1 - means.ifr) / means.ifr))
    pars.ifr[13] <- 9.50
    pars.ifr[14] <- 112
}

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
if(grepl("Anne", scenario.name, fixed = T)){
    ddelay.mean <-15.0
    ddelay.sd <- 12.1
} else {
    ddelay.mean <- 17.8
    ddelay.sd <- 8.9
}

## Reporting delay on the deaths
## First, write down Tom's cdf for the delay distribution function
F <- c(0, 0.0581397, 0.3410914, 0.5708017, 0.7027247, 0.7794353, 0.8213025, 0.8516836, 0.8770857, 0.8940002, .9123456, 0.924536, .9387407, .9544153, .9685447, .9821035, .9868954, .991038, 1)
f <- diff(F)
delay <- 0:(length(f) - 1)
m <- sum(f * delay)
v <- sum(f * delay^2) - m^2
if(grepl("reports", data.desc, fixed = T)){
    rdelay.mean <- m # 3.5
    rdelay.sd <- sqrt(v) # 1.75
} else {
    rdelay.mean <- 0
    rdelay.sd <- 0
}
## Delays in the line listing data
ldelay.mean <- 5.22
ldelay.sd <- 3.59

## Contact model
if(nA == 1){
prior.list <- list(fixed = NULL,
                   tight = c(20.67, 31.00),
                   relax = c(2, 3),
                   unif = c(1, 1),
                   uninf = c(0.5, 0.5))
} else prior.list <- list(relax = c(4, 4))

nprior.name <- names(prior.list)
which.var <- which(sapply(nprior.name, grepl, x = scenario.name, fixed = TRUE))
contact.dist <- rep(c(1, ifelse(is.null(prior.list[which.var]), 1, 3)), nr)
contact.pars <- rep(prior.list[[which.var]], nr)
which.var <- which(sapply(names(int.effect), grepl, x = scenario.name, fixed = TRUE))
stopifnot(length(which.var) == 1)
contact.reduction <- rep(c(1, int.effect[which.var]), nr)
