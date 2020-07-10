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
value.dI <- 0.745258455864861
pars.dI <- c(1.43, 0.549)

## Exponential growth rate
value.egr <- c(0.266386412775556, 0.231040986029484, 0.235589106443936, 0.305950653647048, 0.282526430208455, 0.216668361008301, 0.222294770391012)
pars.egr <- c(31.36, 224)

## Ascertainment parameters
value.pgp <- 0.1
pars.pgp <- c(2.12, 15.8)

## Infection to fatality ratio
if(nA > 1){
    value.ifr <- c(6.35752195427293e-06, 1.25584130296705e-05, 4.07073854458609e-05, 0.000266842164572681, 0.00418815604437853, 0.0280403506685224, 0.22320535744018)
    var.ifr <- rep(0.000368, nA - 1)
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
value.nu <- c(-16.2719920280182, -13.7528046318741, -15.0958323138269, -17.6712207723268, -16.6386535187948, -14.9096259051349, -15.9902773129261)
pars.nu <- c(-17.5, 1.25)

## GP Overdispersion
value.eta <- 0.0406967;
pars.eta <- c(1.0, 0.2);

## Hosp Overdispersion
value.eta.h <- 0.10878003358794
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
contact.proposal <- c(0.0, 0.003313, 0.000031, 0.000564, 0.0, 0.000714, 0.000077, 0.002036, 0.0, 0.001323, 0.000022, 0.000406, 0.0, 0.008052, 0.000020, 0.000216, 0.0, 0.004521, 0.000023, 0.000421, 0.0, 0.001625, 0.000035, 0.000356, 0.0, 0.005989, 0.000097, 0.000872)
## contact.proposal <- rep(c(0, rep(0.0001, nm)), nr)
contact.reduction <- rnorm(length(contact.proposal))
contact.reduction[1 + ((nm+1)*(0:(nr-1)))] <- 0
contact.reduction[contact.reduction != 0] <- c(-0.116541775881024, -0.408980375330183, -0.0156679133877321, -0.0748016516321573, -0.715465768996106, -0.818977354480832, -0.107390822382017, -0.453187349344851, -0.101088273069055, -0.302898386182425, -0.456799393816269, 0.10949140390645, -0.215978158299037, -0.133595830062586, 0.231109645064602, -0.222289432134629, -0.438074600061848, 0.111207936928841, -0.197752353677324, -0.057059948082588, 0.351057759491446)
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
beta.rw.vals <- c(0, -0.0642613930028543, 0.137470227612429, 0.18203661462353, 0.208765769249931, -0.285692841926755, -0.0299902301756172, -0.137239973812927, -0.128225668742848, -0.248624010808639, 0, 0, 0,
                  0, 0.0730655247771195, 0.365275569914078, 0.170953206999095, 0.34635260810452, -0.0802801069213467, -0.195542202819357, -0.0250765553923469, -0.114324681894421, 0.0810947664489636, 0, 0, 0,
                  0, 0.0772995759950353, 0.247531598218536, 0.122210787051954, 0.0444304457953928, -0.00327043114167715, -0.134942895712959, -0.17687660240596, 0.272007780756127, -0.0372792500753984, 0, 0, 0,
                  0, -0.0454182177307798, -0.0346139304577103, 0.360265801369638, 0.0793595872920018, -0.190393947518045, -0.159883554584004, -0.0696252523038599, -0.0549656033022808, -0.0418042501378981, 0, 0, 0,
                  0, -0.33687641540695, -0.0933680963165255, 0.325000771601505, 0.058438676744763, 0.0368240314059391, 0.0229094682674597, 0.0534805446221636, -0.16402496895299, 0.122893481697912, 0, 0, 0,
                  0, 0.212609309121448, 0.041600797861119, 0.115181751556874, 0.0253143971769989, 0.263657560485213, -0.204017906056233, -0.0736330304028773, -0.136342952530373, 0.123952747316385, 0, 0, 0,
                  0, -0.12394074164419, -0.0943171828153658, 0.0563763749693666, -0.0187639431158053, 0.331927977648457, -0.26183947078463, 0.1536281048271, -0.0961824810293474, 0.112915829534447, 0, 0, 0)
beta.update <- TRUE
beta.rw.flag <- TRUE
## beta.rw.props <- rep(c(0, rep(0.0005, nbetas - 1)), nr)
beta.rw.props <- c(0.0, 0.000026, 0.000050, 0.000041, 0.000104, 0.000272, 0.000507, 0.001690, 0.004086, 0.013020, 0.01, 0.01, 0.01,
                   0.0, 0.000064, 0.000057, 0.000147, 0.000216, 0.000398, 0.001232, 0.002377, 0.006850, 0.008003, 0.01, 0.01, 0.01,
                   0.0, 0.000010, 0.000023, 0.000038, 0.000056, 0.000201, 0.000328, 0.000916, 0.002964, 0.009492, 0.01, 0.01, 0.01,
                   0.0, 0.000013, 0.000021, 0.000046, 0.000072, 0.000160, 0.000461, 0.001537, 0.003515, 0.009118, 0.01, 0.01, 0.01,
                   0.0, 0.000013, 0.000019, 0.000032, 0.000098, 0.000124, 0.000218, 0.000880, 0.002066, 0.005020, 0.01, 0.01, 0.01,
                   0.0, 0.000020, 0.000028, 0.000038, 0.000122, 0.000222, 0.000458, 0.000816, 0.004169, 0.010337, 0.01, 0.01, 0.01,
                   0.0, 0.000075, 0.000084, 0.000121, 0.000243, 0.000402, 0.000800, 0.003132, 0.007688, 0.010336, 0.01, 0.01, 0.01)
beta.design <- matrix(1, nbetas, nbetas)
for(i in 1:(nbetas-1))
    for(j in (i+1):nbetas)
        beta.design[i,j] <- 0
beta.design <- as.matrix(bdiag(lapply(1:nr, function(x) beta.design)))
write_tsv(as.data.frame(beta.design), file.path(out.dir, "beta.design.txt"), col_names = FALSE)
beta.rw.sd.pars <- c(1, 100)
beta.rw.sd <- 0.142453430468293

## Serological test sensitivity and specificity
## sero.sens <- 71.5 / 101
## sero.spec <- 777.5 / 787
sero.sens <- 0.707875480848508
sero.spec <- 0.965012479451016
ssens.prior.dist <- ifelse(grepl("var", scenario.name), 3, 1)
## ssens.prior.pars <- c(137.5, 36.5) ## Change the .Rmd file to allow for stochasticity in the sensitivity/specificity
## Default is based on testing intervals 21-27 days, alternative is based on all testing intervals >21 days.
if (grepl("altSens", scenario.name)) {
	ssens.prior.pars <- c(132.5, 28.5)
} else {
	ssens.prior.pars <- c(23.5, 9.5)
}

sspec.prior.dist <- ifelse(grepl("var", scenario.name), 3, 1)
## sspec.prior.pars <- c(699.5, 8.5)
sspec.prior.pars <- c(569.5, 5.5)
ssens.prop <- 0.1
sspec.prop <- 0.077976
