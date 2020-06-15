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
value.dI <- 0.2568197
pars.dI <- c(1.43, 0.549)

## Exponential growth rate
value.egr <- c(0.264283434928516, 0.231964603189632, 0.236287108646892, 0.303878783994188, 0.238121555851514, 0.232136713695436, 0.213344218931685)
pars.egr <- c(31.36, 224)

## Ascertainment parameters
value.pgp <- 0.1
pars.pgp <- c(2.12, 15.8)

## Infection to fatality ratio
if(nA > 1){
    value.ifr <- c(9.75325253605021e-06, 8.75691123158215e-06, 5.39553428326396e-05, 0.000283102986823073, 0.00456517734001891, 0.0296903848028439, 0.154823958035128)
    var.ifr <- rep(0.000365, nA - 1)
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
value.nu <- c(-16.2400138085247, -13.8438735027695, -15.1629831370676, -17.6392471449826, -15.5294374766338, -15.5479280063189, -15.6911571077026)
pars.nu <- c(-17.5, 1.25)

## GP Overdispersion
value.eta <- 0.0406967;
pars.eta <- c(1.0, 0.2);

## Hosp Overdispersion
value.eta.h <- 0.0406967;
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
contact.proposal <- c(0.0, 0.003153, 0.000043, 0.000599, 0.0, 0.000725, 0.000081, 0.001357, 0.0, 0.001493, 0.000026, 0.000261, 0.0, 0.005586, 0.000026, 0.000232, 0.0, 0.004501, 0.000038, 0.000361, 0.0, 0.003520, 0.000026, 0.000354, 0.0, 0.008178, 0.000109, 0.000886)
## contact.proposal <- rep(c(0, rep(0.0001, nm)), nr)
contact.reduction <- rnorm(length(contact.proposal))
contact.reduction[1 + ((nm+1)*(0:(nr-1)))] <- 0
contact.reduction[contact.reduction != 0] <- c(0.210743353677169, -0.325600008144382, 0.473641296269596, 0.396067807844787, -0.712656869139014, -0.394578320242017, 0.15871939824711, -0.298845908557913, 0.546720762440809, 0.0755446120186602, -0.443297508043497, 0.529412723973099, 0.186299939191058, 0.105031817899126, 0.890146864093364, 0.360029530469104, -0.241520652674724, 0.519506891518561, 0.217344189565156, -0.244644487776478, 0.674268621996671)
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
beta.rw.vals <- c(0, -0.152872393737844, 0.197563729934305, 0.334249250865956, -0.0113930340270764, -0.219236314351579, -0.0416932600727287, 0.000531869742348525, 0.0105908128368036, -0.266482286358495, 0, 0.309722784823974, 0.115357295904507, 0.271383321778447, 0.202353477819453, -0.0453584531969326, -0.0773679335946866, 0.173576246204466, 0.0369886543926375, -0.150171389998007, 0, -0.0880096293537438, 0.207429073845858, 0.213588075317703, 0.120197959026136, -0.0678740267942056, -0.121116840018356, -0.00107516915066937, -0.0721185472094092, 0.0653581923735481, 0, 0.0205020541222985, 0.19587316421062, -0.0531717074075164, 0.217465048595645, -0.151944848274846, -0.135496907037351, -0.02284149920891, -0.0020222405646863, -0.21758391298888, 0, -0.188657817126277, -0.291998715622014, 0.279847774045288, 0.155042313997516, 0.0589320634102544, 0.00891460693724563, -0.186493324838978, 0.132666421884143, -0.0178016133164924, 0, 0.0430909859579732, 0.0937771039983638, 0.0225545925249073, 0.174664947484495, 0.0475067447092583, -0.138627981711317, -0.0183209936300398, 0.00949582395629406, -0.104419061109281, 0, 0.141858991340665, 0.0086596842611849, -0.0231104380077788, 0.129476257601218, 0.0148621449435352, -0.0291477915937377, 0.117116984755575, -0.114844867795275, 0.0122794945930524)
beta.update <- TRUE
beta.rw.flag <- TRUE
## beta.rw.props <- rep(c(0, rep(0.0005, nbetas - 1)), nr)
beta.rw.props <- c(0.000000, 0.000026, 0.000037, 0.000070, 0.000219, 0.000364, 0.001077, 0.002590, 0.005693, 0.012144, 0.000000, 0.000085, 0.000085, 0.000154, 0.000206, 0.000520, 0.001349, 0.003279, 0.006718, 0.016163, 0.000000, 0.000016, 0.000024, 0.000025, 0.000092, 0.000106, 0.000459, 0.001433, 0.007693, 0.018696, 0.000000, 0.000013, 0.000018, 0.000039, 0.000077, 0.000195, 0.000413, 0.001221, 0.004497, 0.017602, 0.000000, 0.000024, 0.000039, 0.000028, 0.000115, 0.000156, 0.000482, 0.001104, 0.003553, 0.012144, 0.000000, 0.000016, 0.000036, 0.000079, 0.000116, 0.000246, 0.000552, 0.001946, 0.003535, 0.016327, 0.000000, 0.000056, 0.000145, 0.000172, 0.000224, 0.000692, 0.000894, 0.003379, 0.007316, 0.017602)
beta.design <- matrix(1, nbetas, nbetas)
for(i in 1:(nbetas-1))
    for(j in (i+1):nbetas)
        beta.design[i,j] <- 0
beta.design <- as.matrix(bdiag(lapply(1:nr, function(x) beta.design)))
write_tsv(as.data.frame(beta.design), file.path(out.dir, "beta.design.txt"), col_names = FALSE)
beta.rw.sd.pars <- c(1, 100)
beta.rw.sd <- 0.138219970744751

## Serological test sensitivity and specificity
## sero.sens <- 71.5 / 101
## sero.spec <- 777.5 / 787
sero.sens <- 0.764009681048664
sero.spec <- 0.965192460834826
ssens.prior.dist <- ifelse(grepl("var", scenario.name), 3, 1)
## ssens.prior.pars <- c(137.5, 36.5) ## Change the .Rmd file to allow for stochasticity in the sensitivity/specificity
ssens.prior.pars <- c(23.5, 9.5)
sspec.prior.dist <- ifelse(grepl("var", scenario.name), 3, 1)
## sspec.prior.pars <- c(699.5, 8.5)
sspec.prior.pars <- c(569.5, 5.5)
ssens.prop <- 0.001
sspec.prop <- 0.08
