## Incubation period - best working estimate - mean 5.2 (4.1-7.0)
## Use these as simulation parameters for the latent period

int.effect <- c(0.00, 0.521, 0.24, 1.0)
names(int.effect) <- c("total", "variable", "lshtm", "nothing")

## shape.dL <- 35.1
## rate.dL <- 6.76  ## These values give the desired mean with a variance of 0.768
value.dl <- 1

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
ddelay.mean <-15.0
ddelay.sd <- 12.1
if(grepl("newOtoD", scenario.name, fixed = TRUE)){
    ddelay.mean <- 9.3
    ddelay.sd <- 9.7
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
int.effect <- 0.0521
nm <- max(mult.order)
sd <- sqrt(log(5) - 2*log(2))
design.flag <- grepl("2mrw", scenario.name, fixed = TRUE)
rw.flag <- grepl("allmrw", scenario.name, fixed = TRUE)
prior.list <- list(lock = c(log(2) - 0.5*log(5), sd),
                   effect = c(0, sqrt(1/20) * sd),
                   increments = c(0, sqrt(1/20) * sd)
                   )
if(design.flag) prior.list$lock[2] * sqrt(19/20)
## prior.list <- list(relax = c(4, 4))
contact.dist <- rep(c(1, rep(4, nm)), nr)
## contact.pars <- rep(prior.list[[1]], nr)
contact.pars <- array(0, dim = c(2, nm, nr))
contact.pars[, 1, ] <- prior.list$lock
if(nm > 1){
    for(j in 2:nm)
        contact.pars[, j, ] <- prior.list$increments
}
contact.proposal <- rep(c(0, rep(0.001, nm)), nr)
contact.reduction <- rep(c(0, jitter(rep(int.effect, nm))), nr)
contact.link <- as.integer(any(contact.dist == 4))
require(Matrix)
if(design.flag){
    sub.design <- matrix(c(1, 1, -1, 1), 2, 2)
}
if(rw.flag){
    sub.design <- matrix(1, nm, nm)
    for(i in 1:(nm-1))
        for(j in (i+1):nm)
            sub.design[i,j] <- 0
}
if(design.flag | rw.flag){
    m.design <- lapply(1:nr, function(x) bdiag(1, sub.design))
    m.design <- bdiag(m.design)
    write_tsv(as.data.frame(as.matrix(m.design)), file.path(out.dir, "m.design.txt"), col_names = FALSE)
}

## Serological test sensitivity and specificity
sero.sens <- 71.5 / 101
sero.spec <- 777.5 / 787

ssens.prior.dist <- ifelse(grepl("var", scenario.name), 3, 1)
ssens.prior.pars <- c(71.5, 29.5) ## Change the .Rmd file to allow for stochasticity in the sensitivity/specificity
sspec.prior.dist <- ifelse(grepl("var", scenario.name), 3, 1)
sspec.prior.pars <- c(777.5, 9.5)

ssens.prop <- 0.001
sspec.prop <- 0.001
