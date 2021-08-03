
add.extra.vals.per.region <- function(vec, val, num) {
  mat <- matrix(vec, ncol = nr)
  rows.to.add <- num - length(vec) / nr
  if(rows.to.add > 0){
    mat.new <- matrix(val, nrow = rows.to.add, ncol = nr)
    } else mat.new <- NULL
  return(rbind(mat, mat.new))
}

logit <- function(p) log(p/(1-p))
expit <- function(x) exp(x)/(1+exp(x))

## Incubation period - best working estimate - mean 5.2 (4.1-7.0)
## Use these as simulation parameters for the latent period

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

## Duration of test positivity
value.r1 <- 5.5
prior.r1 <- ifelse(prev.flag, 2, 1)
if(prev.prior == "tight") pars.r1 <- c(550000,100000)
if(prev.prior == "relax") pars.r1 <- c(5.5, 1)
if(prev.prior == "long_positive") pars.r1 <- c(11.7, 0.903)
if(prev.prior == "Cevik") pars.r1 <- c(32.2, 2.60)

vac.effec.bp <- 1:(length(vac.dates) - 1)

## Efficacy against disease from one vaccine dose
if(efficacies == "Nick"){
    value.vac.alpha1 <- c(0.80, 0.50)
} else if(efficacies == "Jamie"){
    value.vac.alpha1 <- rep(3/7,2)
} else if(efficacies == "PHE"){
    value.vac.alpha1 <- c(2/5, 5/14, 31/46, 31/46)
} else {
    value.vac.alpha1 <- c(0.88, 0.70) ## efficacy against disease of Pfizer and AZ vaccines respectively.
}
prior.vac.alpha1 <- rep(1, length(value.vac.alpha1)) ## ifelse(vacc.flag, 3, 1)
prior.alpha1 <- max(prior.vac.alpha1)
if(vacc.flag & (prior.alpha1 > 1)) pars.alpha1 <- c(4, 1)
if(efficacies == "PHE"){
    delta <- as_tibble(expand.grid(age.labs, vac.dates, regions)) %>% mutate(delta = Var2 > delta.date) %>% pull(delta)
    v1.design <- cbind(v1.design, 0, 0)
    v1.design[delta, 3:4] <- v1.design[delta, 1:2]
    v1.design[delta, 1:2] <- 0
}
write_tsv(as.data.frame(v1.design), file.path(out.dir, "vac.alpha1.design.txt"), col_names = FALSE)
vacc.alpha.bps <- TRUE

## Efficacy against disease from second/nth vaccine dose
if(efficacies == "Nick"){
    value.vac.alpha2 <- c(0.95, 0.70)
} else if(efficacies == "Jamie"){
    value.vac.alpha2 <- c(2/3,6/7)
} else if(efficacies == "PHE"){
    value.vac.alpha2 <- c(17/20, 51/57, 17/20, 17/20)
} else {
    value.vac.alpha2 <- c(0.94, 0.82) ## efficacy against disease of Pfizer and AZ vaccines respectively.
}
prior.vac.alpha2 <- rep(1, length(value.vac.alpha2)) ## ifelse(vacc.flag, 3, 1)
prior.alpha2 <- max(prior.vac.alpha2)
if(vacc.flag & (prior.alpha2 > 1)) pars.alpha2 <- c(4, 1)
if(efficacies == "PHE"){
    vn.design <- cbind(vn.design, 0, 0)
    vn.design[delta, 3:4] <- vn.design[delta, 1:2]
    vn.design[delta, 1:2] <- 0
}
write_tsv(as.data.frame(vn.design), file.path(out.dir, "vac.alphan.design.txt"), col_names = FALSE)

## Efficacy against disease from one vaccine dose
if(efficacies == "Nick"){
    value.vac.pi1 <- 0.48 ## This will need to change when I can figure out how(!)
} else if(efficacies == "Jamie"){
    value.vac.pi1 <- c(0.65, 0.65)
} else if(efficacies == "PHE"){
    value.vac.pi1 <- c(0.625, 0.65, 0.31, 0.31)
} else {
    value.vac.pi1 <- 0.48
}
prior.vac.pi1 <- rep(1, length(value.vac.pi1)) ## ifelse(vacc.flag, 3, 1)
prior.pi1 <- max(prior.vac.pi1)
if(vacc.flag & (prior.pi1 > 1)) pars.pi1 <- c(4, 1)
vacc.pi.bps <- efficacies %in% c("Jamie", "PHE")
if(vacc.pi.bps)
    write_tsv(as.data.frame(v1.design), file.path(out.dir, "vac.pi1.design.txt"), col_names = FALSE)

## Efficacy against disease from two vaccine doses
if(efficacies == "Nick"){
    value.vac.pi2 <- 0.6 ## This will need to change when I can figure out how(!)
} else if(efficacies == "Jamie"){
    value.vac.pi2 <- c(0.85, 0.65)
} else if(efficacies == "PHE"){
    value.vac.pi2 <- c(0.8, 0.715, 0.8, 0.8)
} else {
    value.vac.pi2 <- 0.6
}
prior.vac.pi2 <- rep(1, length(value.vac.pi2)) ## ifelse(vacc.flag, 3, 1)
prior.pi2 <- max(prior.vac.pi2)
if(vacc.flag & (prior.pi2 > 1)) pars.pi2 <- c(4, 1)
if(vacc.pi.bps)
    write_tsv(as.data.frame(vn.design), file.path(out.dir, "vac.pin.design.txt"), col_names = FALSE)

## Exponential growth rate
value.egr <- c(0.281224110810985, 0.246300679874443, 0.230259384150778, 0.307383663711624, 0.249492140587071, 0.224509782739688, 0.234528728809235, 0.2, 0.2)[1:nr]
pars.egr <- c(31.36, 224)

## Ascertainment parameters
if(!gp.flag | nA == 1){
    value.pgp <- 0.1
    pars.pgp <- c(2.12, 15.8)
} else if (gp.flag){
    abreaks.icr <- 3:7
    if(grepl("base_varSens", scenario.name)){
        ## For each region and age, get the number of cases immediately prior to the inclusion of the Pillar 2 data.
        ll.prior.days <- ll.start.date - days(1:7)
        ll.prior <- all_dat %>%
            filter(Date %in% ll.prior.days) %>%
            group_modify(~ {.x %>% mutate(Age = ifelse(Age.Grp %in% c("<1yr", "1-4", "<1yr,1-4", "5-14"), "<15yr", Age.Grp))}) %>%
            group_by(Region, Age) %>%
            summarise(npos = sum(npos) / 7,
                      nbar = ifelse("nbar" %in% names(all_dat), sum(nbar) / 7, 1)
                      )
        ## Find an estimate of incidence at this time.
        idx <- which(lubridate::as_date(as.integer(dimnames(outpp$infections)$date)) == min(ll.prior.days))
        infections <- outpp$infections[, idx, , ] ## get infections on the first day of the aggregation
        infections <- apply(infections,
                            c("iteration", "region"),
                            function(x) c(sum(x[1:3]), x[-(1:3)])) ## group first three age groups together
        names(dimnames(infections))[1] <- "age"
        dimnames(infections)$age[1] <- "<15yr"
        med.infec <- apply(infections, c("age", "region"), median) %>%
            as.data.frame() %>%
            rownames_to_column(var = "Age") %>%
            pivot_longer(-Age, names_to = "Region") %>%
            inner_join(ll.prior) %>%
            mutate(p = npos / value) %>%
            arrange(Region, Age)
        ## Weighted regression
        ex6 <- suppressWarnings(glm(p~Region + Age + log(nbar),
                                    weights = value,
                                    family = binomial,
                                    data = med.infec,
                                    contrasts = list(Age = "contr.sum", Region = "contr.sum")))
        ## Use mean and covariance of estimates to formulate a prior distribution.
        if(case.positivity){
            value.pgp <- jitter(ex6$coefficients)
            pars.pgp <- as.vector(t(cbind(ex6$coefficients, round(vcov(ex6), digits = 5))))
            ll.posterior.days <- start.date - 1 + (start.gp:end.gp)
            ll.posterior <- all_dat %>%
                filter(Date %in% ll.posterior.days) %>%
                group_modify(~ {.x %>% mutate(Age.Grp = ifelse(Age.Grp %in% c("<1yr", "1-4", "<1yr,1-4", "5-14"), "<15yr", Age.Grp))}) %>%
                group_by(Region, Date, Age.Grp) %>%
                summarise(n = sum(n),
                          pop = sum(pop)) %>%
                mutate(nbar = n / pop)
            ex6 <- glm(n ~ Region + Age.Grp + log(nbar), data = ll.posterior, contrasts = list(Age.Grp = "contr.sum", Region = "contr.sum"))
            mex6 <- model.matrix(ex6)
        } else {
            nl <- length(ex6$coefficients)
            value.pgp <- jitter(ex6$coefficients[-nl])
            pars.pgp <- as.vector(t(cbind(ex6$coefficients[-nl], round(vcov(ex6)[-nl,-nl], digits = 5))))
            mex6 <- model.matrix(ex6)[,-nl]
        }
    } else {
        reg.mod.loc <- file.path(dirname(proj.dir), "Pillar2", "ascertatinment_regs.RData")
        load(reg.mod.loc)
        idx <- which(is.na(rfull.noinfec$coefficients$x$data[, 2]))
        value.pgp <- rfull.noinfec$coefficients$x$data[-idx, 2]
        pars.pgp <- as.vector(t(cbind(rep(0, times = length(value.pgp)), rep(5, times = length(value.pgp)))))
        mex6 <- rfull.noinfec$design[,-idx]
    }
    write_tsv(as.data.frame(mex6), file.path(out.dir, "icr.design.txt"), col_names = FALSE)
    if(exists("infections"))
        rm(infections)
    if(pgp.prior.diffuse)
        pars.pgp <- as.vector(t(cbind(rep(0, times = length(ex6$coefficients)), rep(2.5, times = length(ex6$coefficients)))))
}

## Infection to fatality ratio
if(single.ifr){
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
} else { ## IFR changes over time, presume logistically...
    if(nA > 1){
        value.ifr <- logit(c(5.77188478860128e-06, 9.57867182705255e-06, 4.5278018816958e-05, 0.000323870211139221, 0.00471791669192503, 0.0316645720110774, 0.202480672513791))
        var.ifr <- rep(0.000360, length(value.ifr))
    } else {
        value.ifr <- logit(c(0.007, 1))
        var.ifr <- rep(0.005, 2)
    }
    if(nA == 1){
        pars.ifr <- c(0.00705, 0.00304)
    } else {
        means.ifr <- c(1.61e-5, 4.28e-5, 1.89e-4, 9.02e-4, 8.20e-3, 3.11e-2, 6.04e-2)
        pars.ifr <- as.vector(rbind(rep(1, length(means.ifr)), (1 - means.ifr) / means.ifr))
        pars.ifr[13] <- 9.50
        pars.ifr[14] <- 112
        ifr <- rbeta(7000000, shape1 = pars.ifr[seq(1, length(pars.ifr), by = 2)], shape2 = pars.ifr[seq(2, length(pars.ifr), by = 2)])
        ifr <- matrix(ifr, nrow = 1000000, ncol = 7, byrow = TRUE)
        pars.ifr <- as.vector(rbind(apply(logit(ifr), 2, mean), apply(logit(ifr), 2, sd)));rm(ifr) ## base parameters - transformed from the informative beta distributions
        if(NHS28.alt.ifr.prior){
            mean.indices <- seq(1, length(pars.ifr), by = 2)
            pars.ifr[mean.indices] <- pars.ifr[mean.indices] - (0.8 * exp(pars.ifr[mean.indices])) + log(0.8)
        }
        ## gradients from Brian's Co-CIN analysis
        grad.samp <- cbind(
        (logit(rbeta(1000000, 23/4.625,977/4.625))-logit(rbeta(1000000, 13.5, 39*13.5))) / 30, ## 0-44
        (logit(rbeta(1000000, 76/2.875, 924/2.875))-logit(rbeta(1000000, 23*5.25, 175*5.25))) / 30, ## 45-64
        (logit(rbeta(1000000, 195/2.9375, 805/2.9375))-logit(rbeta(1000000, 281/1.240235,719/1.240235))) / 30, ## 65-74
        (logit(rbeta(1000000, 321/1.375, 679/1.375))-logit(rbeta(1000000, 382*2.625, 618*2.625))) / 30 ## 75+
        )
        ## expand pars.ifr according to the model
        pars.ifr <- c(pars.ifr, as.vector(apply(grad.samp, 2, function(x) c(mean(x), sd(x)))))
        if(bp.flag <- !is.na(str_extract(ifr.mod, "[0-9]bp"))){
            pars.ifr <- c(pars.ifr, as.vector(apply(grad.samp, 2, function(x) c(-mean(x), sd(x)))))
            if((num.bp <- as.numeric(str_extract(ifr.mod, "[0-9]"))) > 2){
                for(i in num.bp:3)
                    pars.ifr <- c(pars.ifr, as.vector(apply(grad.samp, 2, function(x) c(0, sd(x)))))
            }
        } else if(ifr.mod == "lin.bp"){
            pars.ifr <- c(pars.ifr, as.vector(apply(grad.samp, 2, function(x) 0.3 * c(-mean(x), sd(x)))))
        }; rm(grad.samp)
        value.ifr <- c(value.ifr, rep(0, (length(pars.ifr) / 2) - length(value.ifr)))
        var.ifr <- c(var.ifr, rep(0.00036, length(value.ifr) - length(var.ifr)))
        tbreaks.ifr <- ((ymd("20200531") - start.date):(ymd("20200629") - start.date)) + 1 ## breakpoints from 31st May to 29th June
        ## Put together the model matrix.
        times <- c(tbreaks.ifr, max(tbreaks.ifr) + 1) - min(tbreaks.ifr)
        ages <- 1:(nA - 1)
        if(bp.flag){ ## Use era variable to introduce a second round of breakpoints
            TA <- expand.grid(ages, times, era <- 0:(num.bp - 1)); colnames(TA) <- c("age", "time", "era")
        } else if(ifr.mod == "lin.bp"){
            TA.sub <- expand.grid(ages, times, times.tmp <- 0); colnames(TA.sub) <- c("age", "time1", "time2")
            TA.sub2 <- expand.grid(ages, times.tmp <- max(times), times2 <- 1:100); colnames(TA.sub2) <- c("age", "time1", "time2")
            TA <- rbind(TA.sub, TA.sub2)
        } else {TA <- expand.grid(ages, times);colnames(TA) <- c("age", "time")}
        TA$age.grad <- TA$age
        TA$age.grad[TA$age <= 4] <- 4
        TA$age.grad <- factor(TA$age.grad);TA$age <- factor(TA$age)
        if(bp.flag){ ## Expand the actual breakpoints and tweak the design
            if(!exists("tbreaks.interval")) tbreaks.interval <- min(tbreaks.ifr)
            tbreaks2 <- round(tbreaks.ifr + (rep(1:(num.bp-1), each=length(tbreaks.ifr))*tbreaks.interval) - 1)
            tbreaks.ifr <- c(tbreaks.ifr, tbreaks2)
            reg.form <- "y ~ 0 + age"  ## + age.grad:full.era + age.grad:time:era"
            for(per in 1:num.bp){

                TAcolname <- paste0("time", per)
                TA <- bind_rows(TA %>% filter(era < (per - 1)) %>% mutate(!!TAcolname := 0),
                                TA %>% filter(era == (per - 1)) %>% mutate(!!TAcolname := time),
                                TA %>% filter(era > (per - 1)) %>% mutate(!!TAcolname := max(TA$time))
                                )
                reg.form <- paste0(reg.form,  " + age.grad:", TAcolname)
                if(per > 1){ ## Remove some redundant rows
                    TAprevcolname <- paste0("time", per - 1)
                    TA <- TA %>% filter(!(!!sym(TAprevcolname) == max(TA$time) & !!sym(TAcolname) == 0 & time == 0))
                }
            }

        } else if(ifr.mod == "lin.bp"){
            tbreaks2 <- ymd(date.data) - start.date - (99:0)
            tbreaks.ifr <- c(tbreaks.ifr, tbreaks2)
            reg.form <- "y ~ 0 + age + age.grad:time1 + age.grad:time2"
        } else reg.form <- "y ~ 0 + age + age.grad:time"
        TA$y <- rnorm(nrow(TA))
        lm.TA2 <- lm(as.formula(reg.form), data = TA)
        write_tsv(as.data.frame(model.matrix(lm.TA2)), file.path(out.dir, "ifr.design.txt"), col_names = FALSE)
    }
}

## Day of the week effects in the reporting of `deaths'.
if(gp.flag){
    require(Matrix)
    bank.holiday.days <- ymd(c("20200101", "20200410", "20200413", "20200508", "20200525", "20200831", "20201225", "20201228"))
    ll.days <- ll.start.date + days(0:(end.gp - start.gp))
    DAYS <- wday(ll.days, label = TRUE)
    DAYS[ll.days %in% bank.holiday.days] <- "Sun"
    lm.mat <- model.matrix(~ DAYS, contrasts = list(DAYS = "contr.sum"))[, -1] %>%
        as.data.frame() %>%
        write_tsv(file.path(out.dir, "d_o_w_design_file.txt"), col_names = FALSE)
    if(grepl("base_varSens", scenario.name)){
        pars.dow <- c(0.9840207, 0.9915403, 1.0024927, 0.98667, -0.6, -2.1)
        prop.dow <- rep(2.2e-04, 6)
    } else {
        pars.dow <- rep(0, 6)
        prop.dow <- rep(0, 6)
    }
}

## Initial seeding
value.nu <- c(-16.7064769395683, -14.3327333035582, -15.0542472007424, -17.785749594464, -15.8038617705659, -15.3720269985272, -16.3667281951197, -16, -16)[1:nr]
if(gp.flag){
    if(exists("med.infec")){
        value.nu <- value.nu + logit(med.infec %>% filter(Age == "<15yr") %>% pull(p))
    } else {
        value.nu <- value.nu + logit(0.02)
    }
}
pars.nu <- c(-17.5, 1.25)

## GP Overdispersion
value.eta <- 0.406967;
pars.eta <- c(1.0, 0.02);

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
ldelay.mean <- 3
ldelay.sd <- 3

## ## Contact model ## ## ##
nm <- max(unlist(mult.mat))
sd <- sqrt(log(5) - 2*log(2))
rw.flag <- FALSE
prior.list <- list(lock = c(log(2) - 0.5*log(5), sd),
                   increments = c(0, 0.1 * sd),
                   ## viner = c(-0.775, sqrt(-0.775 - log(0.44))),
                   viner = c(-0.4325, sqrt(-0.4325 - log(0.64))),
                   ons = c(-0.6554688, sqrt(-0.6554688 - log(0.5)))
                   )
## prior.list <- list(relax = c(4, 4))
contact.dist <- rep(c(1, rep(4, nm)), nr)
## contact.pars <- rep(prior.list[[1]], nr)
contact.pars <- array(0, dim = c(2, nm, nr))
for(i in 1:nm){
    contact.pars[, i, ] <- prior.list$lock
    if((contact.model == 4) & (i %in% c(1, 4)))
        contact.pars[, i, ] <- prior.list[[contact.prior]]
    if((contact.model %in% c(5, 6)) & (i %in% c(1, 5)))
        contact.pars[, i, ] <- prior.list[[contact.prior]]
    if((contact.model %in% c(5, 6)) & (i %in% c(2, 6)))
        contact.pars[, i, ] <- 0.5 * (prior.list[[contact.prior]] + prior.list$lock)
}
## if(nm > 1){
##     for(j in 2:nm)
##         contact.pars[, j, ] <- prior.list$increments
## }
contact.proposal <- c(0.001, 0.003462, 0.000018, 0.000463,
                      0.001, 0.000769, 0.000037, 0.001446,
                      0.001, 0.001795, 0.000009, 0.000360,
                      0.001, 0.005577, 0.000010, 0.000311,
                      0.001, 0.003568, 0.000015, 0.000486,
                      0.001, 0.002389, 0.000012, 0.000359,
                      0.001, 0.004945, 0.000047, 0.000871)
contact.proposal <- c(contact.proposal, rep(c(0.001, 0.004945, 0.000047, 0.000871), nr - 7))
if(contact.model == 4)
    contact.proposal <- as.vector(t(matrix(contact.proposal, nr, length(contact.proposal) / nr, byrow = TRUE)[, c(1, 2, 2, 3, 4, 4)]))
if(contact.model == 5)
    contact.proposal <- as.vector(t(matrix(contact.proposal, nr, length(contact.proposal) / nr, byrow = TRUE)[, c(1, 2, 2, 2, 3, 4, 4, 4)]))
if(contact.model == 6)
    contact.proposal <- as.vector(t(matrix(contact.proposal, nr, length(contact.proposal) / nr, byrow = TRUE)[, c(1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4)]))
## contact.proposal <- rep(c(0, rep(0.0001, nm)), nr)
contact.reduction <- rep(-0.1, length(contact.proposal))
zero.contact.elements <- 1 + ((nm+1)*(0:(nr-1)))
contact.reduction[zero.contact.elements] <- 0
contact.proposal[zero.contact.elements] <- 0
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
## ## End of contact model ##

## beta.breaks <- cm.breaks[cm.breaks <= (ymd(date.data) - start.date - time.to.last.breakpoint)][-1] ## Josh's version
beta.breaks.full <- cm.breaks[cm.breaks <= (ndays - (7 * nforecast.weeks) - time.to.last.breakpoint)][-1] ## Paul's version
beta.breaks <- beta.breaks.full[rev(seq(length(beta.breaks.full), 1, by = -break.window))]
beta.inds <- beta.breaks.full %in% beta.breaks

nbetas <- length(beta.breaks) + 1
nbetas.full <- length(beta.breaks.full) + 1
beta.rw.vals <- c(
    0, 0.157278692089345, 0.100224366819243, -0.00751722586102982, 0.242453478585049, 0.115532807450752, -0.0739340379686089, -0.253747153954472, -0.154515402774888, -0.0116745847409388, 0.178320412236642, -0.121550440418646, 0.136056150333538, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, 0.0163618520524716, 0.333494806429032, 0.21527031558633, 0.161083514287401, 0.0463132886726999, -0.0332667493328686, 0.0912165703264888, -0.0674068339766244, -0.191029847903661, -0.289385771665435, -0.0712357129236337, 0.17504363429134, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, -0.00579603465388075, 0.417311445648212, 0.126514200748151, -0.0263963964981505, -0.220355017816457, 0.118815574597186, 0.0313276812519405, -0.0517133569835031, -0.0767156700375162, -0.0381946610630543, -0.231271157270191, -0.13303507053379, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, -0.196305061435333, 0.284216630972756, 0.188116721872681, -0.0342462506148033, -0.0573009315823325, -0.0900636672344598, -0.00425473562060216, -0.32540002766214, 0.104877868600528, 0.114177956716661, 0.0821584040758182, -0.417164459028684, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, -0.348925402477912, 0.0679687980318762, 0.213164755977372, 0.074591979115163, 0.0750080453232727, 0.102336034443526, -0.229377011561744, 0.0520233979921866, -0.279980273239356, 0.0347612260908936, -0.114458350262099, 0.0545339652877542, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, 0.0448275470941772, 0.0513151373848244, 0.0120395022862853, 0.0486208080647384, 0.237665958394784, -0.112122908685769, 0.000419907134729215, -0.0739860667978034, -0.143566919550603, -0.182386385950509, 0.250466537490249, -0.0211042287438713, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, 0.0743603245592828, -0.135251090010906, -0.0360794056507664, 0.110415684736955, 0.109741332977249, 0.155427165123845, -0.0848892480165284, -0.100112415417403, -0.351786922834953, -0.239464175187904, 0.186487858627732, -0.121900557631279, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
)[1:(nbetas.full*nr)]
beta.rw.vals <- add.extra.vals.per.region(beta.rw.vals, 0.0, nbetas.full)
if(length(beta.rw.vals) > nbetas*nr)
    beta.rw.vals <- beta.rw.vals[c(1, 1+sort(sample.int(nbetas.full-1, nbetas-1))),]
static.zero.beta.locs <- seq(from = 1, by = nbetas, length = nr)
beta.dist <- rep(4, length(beta.rw.vals))
beta.dist[static.zero.beta.locs] <- 1
beta.update <- TRUE
beta.rw.flag <- TRUE
## beta.rw.props <- rep(c(0, rep(0.0005, nbetas - 1)), nr)
beta.rw.props <- c(
    0.000000, 0.000011, 0.000012, 0.000020, 0.000038, 0.000053, 0.000090, 0.000114, 0.000344, 0.000958, 0.001271, 0.005583, 0.011967, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000019, 0.000047, 0.000067, 0.000065, 0.000125, 0.000161, 0.000290, 0.000677, 0.001589, 0.003719, 0.005959, 0.010296, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000003, 0.000007, 0.000011, 0.000016, 0.000020, 0.000059, 0.000101, 0.000157, 0.000436, 0.001271, 0.004388, 0.008550, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000005, 0.000008, 0.000008, 0.000023, 0.000030, 0.000063, 0.000131, 0.000279, 0.000616, 0.001303, 0.004852, 0.009842, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000006, 0.000008, 0.000012, 0.000023, 0.000027, 0.000056, 0.000119, 0.000272, 0.000384, 0.001462, 0.004730, 0.013297, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000012, 0.000015, 0.000022, 0.000034, 0.000033, 0.000079, 0.000138, 0.000316, 0.000654, 0.001155, 0.003663, 0.012773, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    0.000000, 0.000018, 0.000033, 0.000067, 0.000120, 0.000247, 0.000245, 0.000519, 0.001663, 0.002612, 0.006823, 0.007852, 0.010879, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02
)
beta.rw.props <- add.extra.vals.per.region(beta.rw.props, 0.02, nbetas.full)
if (length(beta.rw.props) > nbetas.full*nr) beta.rw.props <- beta.rw.props[nbetas.full*nr]
beta.design <- matrix(1, nbetas, nbetas)
for(i in 1:(nbetas-1))
    for(j in (i+1):nbetas)
        beta.design[i,j] <- 0
beta.design <- as.matrix(bdiag(lapply(1:nr, function(x) beta.design)))
write_tsv(as.data.frame(beta.design), file.path(out.dir, "beta.design.txt"), col_names = FALSE)
beta.rw.sd.pars <- c(1, sdpar)
beta.rw.sd <- 0.151057317190954

## Serological test sensitivity and specificity
## sero.sens <- 71.5 / 101
## sero.spec <- 777.5 / 787
sero.sens <- 0.707875480848508
sero.spec <- 0.965012479451016
ssens.prior.dist <- ifelse(fix.sero.test.spec.sens, 1, 3)
## ssens.prior.pars <- c(137.5, 36.5) ## Change the .Rmd file to allow for stochasticity in the sensitivity/specificity
## Default is based on testing intervals 21-27 days, alternative is based on all testing intervals >21 days.
## if (grepl("altSens", scenario.name)) {
## 	ssens.prior.pars <- c(142.5, 29.5)
## 	sspec.prior.pars <- c(1110.5, 8.5)
## } else {
	ssens.prior.pars <- c(23.5, 9.5)
	sspec.prior.pars <- c(569.5, 5.5)
## }

sspec.prior.dist <- ssens.prior.dist
## sspec.prior.pars <- c(699.5, 8.5)
ssens.prop <- 0.1
sspec.prop <- 0.077976
if(ssens.prior.dist == 1) sero.sens <-  ssens.prior.pars[1] / (sum(ssens.prior.pars))
if(sspec.prior.dist == 1) sero.spec <-  sspec.prior.pars[1] / (sum(sspec.prior.pars))

if(use.previous.run.for.start) {
    previous.loc <- previous.run.to.use[1]
    source(file.path(proj.dir, "import_pars.R"))
}

source(file.path(proj.dir, "par_check.R"))
