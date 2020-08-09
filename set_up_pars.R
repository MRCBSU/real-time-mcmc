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
## beta.breaks <- c(64, 92, 120)
nbetas <- length(beta.breaks) + 1
beta.rw.vals <- c(
    0, 0.157278692089345, 0.100224366819243, -0.00751722586102982, 0.242453478585049, 0.115532807450752, -0.0739340379686089, -0.253747153954472, -0.154515402774888, -0.0116745847409388, 0.178320412236642, -0.121550440418646, 0.136056150333538, 0.0, 0.0, 0.0, 0.0,
    0, 0.0163618520524716, 0.333494806429032, 0.21527031558633, 0.161083514287401, 0.0463132886726999, -0.0332667493328686, 0.0912165703264888, -0.0674068339766244, -0.191029847903661, -0.289385771665435, -0.0712357129236337, 0.17504363429134, 0.0, 0.0, 0.0, 0.0,
    0, -0.00579603465388075, 0.417311445648212, 0.126514200748151, -0.0263963964981505, -0.220355017816457, 0.118815574597186, 0.0313276812519405, -0.0517133569835031, -0.0767156700375162, -0.0381946610630543, -0.231271157270191, -0.13303507053379, 0.0, 0.0, 0.0, 0.0, 
    0, -0.196305061435333, 0.284216630972756, 0.188116721872681, -0.0342462506148033, -0.0573009315823325, -0.0900636672344598, -0.00425473562060216, -0.32540002766214, 0.104877868600528, 0.114177956716661, 0.0821584040758182, -0.417164459028684, 0.0, 0.0, 0.0, 0.0,
    0, -0.348925402477912, 0.0679687980318762, 0.213164755977372, 0.074591979115163, 0.0750080453232727, 0.102336034443526, -0.229377011561744, 0.0520233979921866, -0.279980273239356, 0.0347612260908936, -0.114458350262099, 0.0545339652877542, 0.0, 0.0, 0.0, 0.0,
    0, 0.0448275470941772, 0.0513151373848244, 0.0120395022862853, 0.0486208080647384, 0.237665958394784, -0.112122908685769, 0.000419907134729215, -0.0739860667978034, -0.143566919550603, -0.182386385950509, 0.250466537490249, -0.0211042287438713, 0.0, 0.0, 0.0, 0.0, 
    0, 0.0743603245592828, -0.135251090010906, -0.0360794056507664, 0.110415684736955, 0.109741332977249, 0.155427165123845, -0.0848892480165284, -0.100112415417403, -0.351786922834953, -0.239464175187904, 0.186487858627732, -0.121900557631279, 0.0, 0.0, 0.0, 0.0
)[1:(nbetas*nr)]
>>>>>>> 38390767965eeab3309fca53af0a491fb1276c83
beta.update <- TRUE
beta.rw.flag <- TRUE
## beta.rw.props <- rep(c(0, rep(0.0005, nbetas - 1)), nr)
beta.rw.props <- c(
    0.000000, 0.000011, 0.000012, 0.000020, 0.000038, 0.000053, 0.000090, 0.000114, 0.000344, 0.000958, 0.001271, 0.005583, 0.011967, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000019, 0.000047, 0.000067, 0.000065, 0.000125, 0.000161, 0.000290, 0.000677, 0.001589, 0.003719, 0.005959, 0.010296, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000003, 0.000007, 0.000011, 0.000016, 0.000020, 0.000059, 0.000101, 0.000157, 0.000436, 0.001271, 0.004388, 0.008550, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000005, 0.000008, 0.000008, 0.000023, 0.000030, 0.000063, 0.000131, 0.000279, 0.000616, 0.001303, 0.004852, 0.009842, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000006, 0.000008, 0.000012, 0.000023, 0.000027, 0.000056, 0.000119, 0.000272, 0.000384, 0.001462, 0.004730, 0.013297, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000012, 0.000015, 0.000022, 0.000034, 0.000033, 0.000079, 0.000138, 0.000316, 0.000654, 0.001155, 0.003663, 0.012773, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000018, 0.000033, 0.000067, 0.000120, 0.000247, 0.000245, 0.000519, 0.001663, 0.002612, 0.006823, 0.007852, 0.010879, 0.02, 0.02, 0.02, 0.02
)[1:(nbetas*nr)]
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
