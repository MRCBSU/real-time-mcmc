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
value.dI <- 0.5033387
pars.dI <- c(1.43, 0.549)

## Exponential growth rate
value.egr <- c(0.281224110810985, 0.246300679874443, 0.230259384150778, 0.307383663711624, 0.249492140587071, 0.224509782739688, 0.234528728809235)[1:nr]
pars.egr <- c(31.36, 224)

## Ascertainment parameters
value.pgp <- 0.1
pars.pgp <- c(2.12, 15.8)

## Infection to fatality ratio
if(nA > 1){
    value.ifr <- c(5.77188478860128e-06, 9.57867182705255e-06, 4.5278018816958e-05, 0.000323870211139221, 0.00471791669192503, 0.0316645720110774, 0.202480672513791)
    var.ifr <- rep(0.000360, nA - 1)
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
value.nu <- c(-16.7064769395683, -14.3327333035582, -15.0542472007424, -17.785749594464, -15.8038617705659, -15.3720269985272, -16.3667281951197)[1:nr]
pars.nu <- c(-17.5, 1.25)

## GP Overdispersion
value.eta <- 0.0406967;
pars.eta <- c(1.0, 0.2);

## Hosp Overdispersion
value.eta.h <- 0.06729983
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
nm <- max(unlist(mult.mat))
sd <- sqrt(log(5) - 2*log(2))
rw.flag <- FALSE
prior.list <- list(lock = c(log(2) - 0.5*log(5), sd),
                   increments = c(0, 0.1 * sd)
                   )
## prior.list <- list(relax = c(4, 4))
contact.dist <- rep(c(1, rep(4, nm)), nr)
## contact.pars <- rep(prior.list[[1]], nr)
contact.pars <- array(0, dim = c(2, nm, nr))
for(i in 1:nm)
contact.pars[, i, ] <- prior.list$lock
## if(nm > 1){
##     for(j in 2:nm)
##         contact.pars[, j, ] <- prior.list$increments
## }
contact.proposal <- c(0.0, 0.003462, 0.000018, 0.000463, 0.0, 0.000769, 0.000037, 0.001446, 0.0, 0.001795, 0.000009, 0.000360, 0.0, 0.005577, 0.000010, 0.000311, 0.0, 0.003568, 0.000015, 0.000486, 0.0, 0.002389, 0.000012, 0.000359, 0.0, 0.004945, 0.000047, 0.000871)[1:length(contact.dist)]
## contact.proposal <- rep(c(0, rep(0.0001, nm)), nr)
contact.reduction <- rnorm(length(contact.proposal))
contact.reduction[1 + ((nm+1)*(0:(nr-1)))] <- 0
if (nr == 7) {
  contact.reduction[contact.reduction != 0] <- c(-0.0341814316217191, -0.603623169409033, 0.0539584598346591, 0.163593712131155, -0.703814518195703, -0.550731233431052, 0.0504542512327084, -0.420094076795744, 0.204225185679588, -0.324308107735709, -0.434259090524466, 0.380906862841831, -0.196723716378284, -0.0485712125740307, 0.667651849590346, 0.0768011654865117, -0.237594065172005, 0.369835713299122, -0.021905891919917, -0.20964475571682, 0.445364908401734)
}
contact.link <- as.integer(any(contact.dist == 4))
require(Matrix)
if(rw.flag){
    sub.design <- matrix(1, nm, nm)
    for(i in 1:(nm-1))
        for(j in (i+1):nm)
            sub.design[i,j] <- 0
}
if(rw.flag){
    m.design <- lapply(1:nr, function(x) bdiag(1, sub.design))
    m.design <- bdiag(m.design)
    m.design <- as.matrix(m.design)
    write_tsv(as.data.frame(m.design), file.path(out.dir, "m.design.txt"), col_names = FALSE)
}

beta.breaks <- cm.breaks[-c(1, length(cm.breaks))]
nbetas <- length(beta.breaks) + 1
beta.rw.vals <- rep(c(0, jitter(rep(0, nbetas - 1), amount = 0.2)), nr)
beta.update <- TRUE
beta.rw.flag <- TRUE
beta.rw.props <- rep(c(0, rep(5e-4, nbetas - 1)), nr)
beta.design <- matrix(1, nbetas, nbetas)
for(i in 1:(nbetas-1))
    for(j in (i+1):nbetas)
        beta.design[i,j] <- 0
beta.design <- as.matrix(bdiag(lapply(1:nr, function(x) beta.design)))
write_tsv(as.data.frame(beta.design), file.path(out.dir, "beta.design.txt"), col_names = FALSE)
beta.rw.sd.pars <- c(1, 100)
beta.rw.sd <- 0.151057317190954

## Serological test sensitivity and specificity
## sero.sens <- 71.5 / 101
## sero.spec <- 777.5 / 787
sero.sens <- 0.707875480848508
sero.spec <- 0.965012479451016
ssens.prior.dist <- ifelse(grepl("altSens", scenario.name) | grepl("varSens", scenario.name), 3, 1)
## ssens.prior.pars <- c(137.5, 36.5) ## Change the .Rmd file to allow for stochasticity in the sensitivity/specificity
## Default is based on testing intervals 21-27 days, alternative is based on all testing intervals >21 days.
if (grepl("altSens", scenario.name)) {
	ssens.prior.pars <- c(142.5, 29.5)
	sspec.prior.pars <- c(1110.5, 8.5)
} else {
	ssens.prior.pars <- c(23.5, 9.5)
	sspec.prior.pars <- c(569.5, 5.5)
}

sspec.prior.dist <- ssens.prior.dist
## sspec.prior.pars <- c(699.5, 8.5)
ssens.prop <- 0.1
sspec.prop <- 0.077976
