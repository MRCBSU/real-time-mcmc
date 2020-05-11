pw.cut.columns <- function(breakpoints, max.days, intervals.per.day, start = 0, end = max.days){

  cutpoints <- c(0, breakpoints, max.days)
  
  as.numeric(cut(start + (1:((end - start) * intervals.per.day) / intervals.per.day), breaks = cutpoints))
}

### Scaling contact matrices by dominant eigenvalues
ngm.matrix <- function(M, popn){
  NGM <- matrix(0, nrow(M), ncol(M))
  for(j in 1:ncol(M))
    NGM[, j] <- M[, j] * popn
  NGM
}
scaled.mixing.matrix <- function(mixing.model, popn){

  M <- mixing.model$M[[1]]

  NGM <- ngm.matrix(M, popn)
  
  NGM.eval <- max(abs(eigen(NGM)$values))
  
  out <- lapply(mixing.model$M, function(mat) mat / NGM.eval)
}
scaled.init.age.distribution <- function(M){
  initial.age.distribution <- eigen(M)$vectors[, 1]
  if(any(zapsmall(Im(initial.age.distribution)) > 0)){
      stop("Non-real initial age distribution")
      return(-1)
  }
  initial.age.distribution <- Re(initial.age.distribution)
  initial.age.distribution / sum(initial.age.distribution)
}

fn.background <- function(bg.model, nages, ndays, nsteps.day, population.size)
  {
    bg.intervals <- pw.cut.columns(bg.model$breakpoints, ndays, nsteps.day)
    bg.GP.delta.t <- exp(bg.model$design.matrix %*% bg.model$basic.rates)
    bg.GP.delta.t <- matrix(bg.GP.delta.t[bg.model$parameterisation.matrix[, bg.intervals]], nages, length(bg.intervals))
    t(bg.GP.delta.t * population.size / (100000 * 365.25 * nsteps.day))
  }
### FUNCTION TO CARRY OUT THE CONVOLUTION OF INFECTIONS OVER THE DELAY DISTRIBUTIONS TO GET PREDICTED GP CONSULTATIONS ###
reporting.model <- function(NNI, q, proportion.symptomatics, reporting.param.list, proportion.model = NULL, background.model = NULL, dow.effects = NULL, delay.distribution = NULL, intervals.per.day = 2, num.days = 200, population.size)
  {
  
    ## NNI should be a (num.days x intervals.per.day) x num.ages matrix
    num.age.groups <- ncol(NNI)
    
    ## delay.distribution should be set a priori if not dependent on any update parameters.
    if(is.null(delay.distribution))
      delay.distribution <- discretised.delay.cdf(reporting.param.list, intervals.per.day)
    ## browser()
    ## a fixed proportion of new infections have onset of symptoms
    NNI.discount <- NNI * proportion.symptomatics
    
    ## GP consultations and proportion per half day
    GP.delta.t <- matrix(0, (intervals.per.day * num.days), num.age.groups)
    
    if(!is.null(proportion.model)){
      p.intervals <- pw.cut.columns(proportion.model$breakpoints, num.days, intervals.per.day)
      
      p <- unlist(eval(parse(text = proportion.model$p.transformation)))
      
      proportion <- t(matrix(p[proportion.model$parameterisation.matrix[, p.intervals]], num.age.groups, num.days * intervals.per.day))
    } else proportion <- q
    
    if(!is.null(background.model))
        bg.GP.delta.t <- fn.background(background.model, num.age.groups, num.days, intervals.per.day, population.size)
    ## ##################################################################### #

    for(a in 1:num.age.groups)
        GP.delta.t[, a] <- conv(NNI.discount[, a], delay.distribution)[1:nrow(GP.delta.t)]

    GP.delta.t <- GP.delta.t * proportion
    
    days <- gl(num.days, intervals.per.day)
    
    GP.day <- by(GP.delta.t, days, function(x) apply(x, 2, sum))
    if(ncol(GP.delta.t) > 1){
      GP.day <- do.call("rbind", GP.day)
    } else {
      GP.day <- matrix(GP.day, length(GP.day), 1)
    }
    
    if(!is.null(background.model))
      {
        bg.GP.day <- by(bg.GP.delta.t, days, function(x) apply(x, 2, sum))
        if(ncol(GP.delta.t) > 1){
          bg.GP.day <- do.call("rbind", bg.GP.day)
        } else {
          bg.GP.day <- matrix(bg.GP.day, length(GP.day), 1)
        }
      } else bg.GP.day = 0

    if(!is.null(dow.effects))
      {
        dow.perday <- as.vector(exp(dow.effects$design.matrix %*% dow.effects$params))
        if(length(dow.perday) == nrow(GP.day)) ## assumes that there are no age effects.
          { ## things are ok
            GP.day <- GP.day * dow.perday
            bg.GP.day <- bg.GP.day * dow.perday
          } else { ## something has gone wrong
            stop("Day of week effects not equal to number of days")
          }
        
      }
    
    list(excess = GP.day, background = bg.GP.day, total = (GP.day + bg.GP.day))
}

## TRANSMISSION KERNEL FOR REED-FROST
calc.p.lambda.row <- function(p.beta, ## a matrix
                              I.1, I.2, ## vectors
                              delta.t = 0.5, relative.susceptibility = 1)
    {
        if(length(p.beta) == 1) p.beta <- matrix(p.beta, 1, 1)
        outvec <- vector("double", nrow(p.beta))
        for(a in 1:nrow(p.beta))
            outvec[a] <- (1 - prod((1 - p.beta[a, ])^(I.1 + (relative.susceptibility * I.2))))
        outvec * delta.t
    }
## TRANSMISSION KERNEL FOR MASS-ACTION
calc.p.lambda.row4 <- function(p.beta,
                               I.1, I.2,
                               delta.t = 0.5, relative.susceptibility = 1)
    {
        if(length(p.beta) == 1) p.beta <- matrix(p.beta, 1, 1)
        outvec <- vector("double", nrow(p.beta))
        for(a in 1:nrow(p.beta))
            outvec[a] <- sum((I.1 + (relative.susceptibility * I.2)) * p.beta[a, ])
        outvec * delta.t
    }


SEEIIR.model.function <- function(phi, R0.amplitude.param = 0, SEIR.param.list, mixing.model, population.size, country.population, intervals.per.day = 2, num.days = 245, tracing.day = 92, start.day = 121, result.by.half.day = F, chasing.threshold = 30){

    attach(phi) ## list of parameters including mm, egr, aip, nu, theta

    region.population <- sum(population.size)
    num.age.groups <- length(population.size)
    num.intervals <- num.days * intervals.per.day
    num.intervals.nontrace <- (tracing.day - 1) * intervals.per.day
    times <- 1:num.intervals
    delta.t <- 1 / intervals.per.day
    
    R.init <- egr * aip * (((egr * SEIR.param.list$latent.period / 2) + 1)^2) / (1 - (1 / (((egr * aip / 2) + 1)^2)))
    
    R0.amplitude <- R0.amplitude.param * R.init / (1 + cos(2 * pi * (start.day - SEIR.param.list$R0.peak.day) / 365.25))
    
    if(R0.amplitude != 0){
        R0 <- R.init + (R0.amplitude * (cos( 2 * pi * ((times * delta.t) + start.day - SEIR.param.list$R0.peak.day) / 365.25) - cos(2 * pi * (start.day - SEIR.param.list$R0.peak.day) / 365.25)))
    } else R0 <- rep(R.init, length(times))

    ## Overall initial number of infectives
    I0.tot <- exp(nu) * aip * region.population / (pgp * R.init)
            
    ## Mixing model
    mm <- c(mm,
            mm[2], ## PHASE 1
            mm[2] + (0.1 * (1 - mm[2])), ## PHASE 2
            mm[2] + (0.3 * (1 - mm[2])), ## PHASE 3
            mm[2] + (0.75 * (1 - mm[2]))) ## PHASE 4
    mixing.model$M <- mixing.model$base.contact.matrices
    for(j in 1:length(mixing.model$base.contact.matrices))
        mixing.model$M[[j]] <- mixing.model$base.contact.matrices[[j]] * matrix(mm[mixing.model$scalants.param[[j]] + 1], n.age, n.age)
    mixing.model$scaled.matrices <- scaled.mixing.matrix(mixing.model, population.size)
    
    ## INITIALISE MODEL
    rho <-  2 * delta.t/ aip
    sigma <- delta.t * 2 / SEIR.param.list$latent.period
    alpha <- exp(egr * delta.t) - 1
    I0 <- I0.tot * mixing.model$initial.age.distribution
    
    I1 <- I0 / (1 + (rho / (alpha + rho)))
    I2 <- I1 * (rho / (alpha + rho))
    E2 <- I1 * ((alpha + rho) / sigma)
    E1 <- E2 * ((alpha + sigma) / sigma)  
    
    NNI <- p.lambda <- S <- I.1 <- I.2 <- E.1 <- E.2 <- R <- Qr <- Qs <- q <- TE.1 <- TE.2 <- TI.1 <- TI.2 <- UE.1 <- UE.2 <- UI.1 <- UI.2 <- matrix(0, num.intervals, num.age.groups)
    
    I.1[1, ] <- I1
    I.2[1, ] <- I2
    E.1[1, ] <- E1
    E.2[1, ] <- E2
    Qr[1, ] <- Qs[1, ] <- q[1, ] <- 0
    TE.1[1, ] <- UE.1[1, ] <- 0
    TE.2[1, ] <- UE.2[1, ] <- 0
    TI.1[1, ] <- UI.1[1, ] <- 0
    TI.2[1, ] <- UI.2[1, ] <- 0
    S[1, ] <- population.size * SEIR.param.list$susceptibility
    R[1, ] <- (population.size * (1 - SEIR.param.list$susceptibility)) - I1 - I2 - E1 - E2
    
    p.beta <- array(0, c(num.intervals, num.age.groups, num.age.groups))
    
    m.max <- min((1:length(mixing.model$intervals))[sapply(1:length(mixing.model$intervals), function(y) num.intervals %in% mixing.model$intervals[[y]])])
    
    for(idx in 1:m.max){
        t.max <- min(num.intervals, mixing.model$intervals[[idx]][length(mixing.model$intervals[[idx]])])
        t.min <- mixing.model$intervals[[idx]][1]
        for(t in t.min:t.max)
            p.beta[t, , ] <- R0[t] * mixing.model$scaled.matrices[[idx]] / aip
    }
    
    p.lambda[1, ] <- calc.p.lambda.row(p.beta[1, , ], I.1[1, ], I.2[1, ], delta.t, SEIR.param.list$relative.susceptibility.I1.to.I2)
    NNI[1, ] <- p.lambda[1, ] * S[1, ]
    
    for(t in 2:num.intervals.nontrace){
        S[t, ] <- S[t - 1, ] * (1 - p.lambda[t -1, ])
        E.1[t, ] <- ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * E.1[t - 1, ]) + (p.lambda[t - 1, ] * S[t - 1, ])
        E.2[t, ] <- ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * E.2[t - 1, ]) + ((2 * delta.t / SEIR.param.list$latent.period) * E.1[t - 1, ])
        I.1[t, ] <- (I.1[t - 1, ] * (1 - (2 * delta.t / aip))) + ((2 * delta.t / SEIR.param.list$latent.period) * E.2[t - 1, ])
        I.2[t, ] <- (I.2[t - 1, ] * (1 - (2 * delta.t / aip))) + (I.1[t - 1, ] * 2 * delta.t / aip)
        R[t, ] <- R[t - 1, ] + (I.2[t - 1, ] * 2 * delta.t / aip)
        NNI[t, ] <- p.lambda[t - 1, ] * S[t - 1, ]
        p.lambda[t, ] <- calc.p.lambda.row(p.beta[t, , ], I.1[t, ], I.2[t, ], delta.t, SEIR.param.list$relative.susceptibility.I1.to.I2)
        Qr[t, ] <- Qs[t, ] <- 0
        TE.1[t, ] <- UE.1[t, ] <- 0
        TE.2[t, ] <- UE.2[t, ] <- 0
        TI.1[t, ] <- UI.1[t, ] <- 0
        TI.2[t, ] <- UI.2[t, ] <- 0
        q[t - 1] <- 0
    }
    for(t in (num.intervals.nontrace+1):num.intervals){
        NNI[t, ] <- p.lambda[t - 1, ] * S[t - 1, ]
        S[t, ] <- S[t - 1, ] - NNI[t, ] - Qs[t - 4, ] + Qs[t - 14, ]
        E.1[t, ] <- (NNI[t, ] * (1 - q[t - 1])) + ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * E.1[t - 1, ])
        new.I1 <- ((2 * delta.t / SEIR.param.list$latent.period) * E.2[t - 1, ])
        E.2[t, ] <- E.2[t - 1, ] - new.I1 + ((2 * delta.t / SEIR.param.list$latent.period) * E.1[t - 1, ])
        new.I2 <- I.1[t - 1, ] * 2 * delta.t / aip
        I.1[t, ] <- I.1[t - 1, ] - new.I2 + new.I1
        I.2[t, ] <- (I.2[t - 1, ] * (1 - (2 * delta.t / aip))) + new.I2
        TE.1[t, ] <- (SEIR.param.list$p.E * q[t - 1, ] * NNI[t, ]) + ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * TE.1[t - 1, ])
        TE.2[t, ] <- ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * TE.2[t - 1, ]) + ((2 * delta.t / SEIR.param.list$latent.period) * TE.1[t - 1, ])
        TI.1[t, ] <- (TI.1[t - 1, ] * (1 - (2 * delta.t / aip))) + ((2 * delta.t / SEIR.param.list$latent.period) * TE.2[t - 1, ])
        TI.2[t, ] <- (TI.2[t - 1, ] * (1 - (2 * delta.t / aip))) + (TI.1[t - 1, ] * 2 * delta.t / aip)
        UE.1[t, ] <- ((1 - SEIR.param.list$p.E) * q[t - 1, ] * NNI[t, ]) + ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * UE.1[t - 1, ])
        new.UI1 <- ((2 * delta.t / SEIR.param.list$latent.period) * UE.2[t - 1, ])
        UE.2[t, ] <- UE.2[t - 1, ] - new.UI1 + ((2 * delta.t / SEIR.param.list$latent.period) * UE.1[t - 1, ])
        new.UI2 <- UI.1[t - 1, ] * 2 * delta.t / aip
        UI.1[t, ] <- UI.1[t - 1, ] - new.UI2 + new.UI1
        UI.2[t, ] <- (UI.2[t - 1, ] * (1 - (2 * delta.t / aip))) + new.UI2
        R[t, ] <- R[t - 1, ] + ((I.2[t - 1, ] + TI.2[t - 1, ] + UI.2[t - 1, ]) * 2 * delta.t / aip) - Qr[t - 4, ] + Qr[t - 14, ]
        ## Use simulation to get a proportion of time that the I states will transmit that the UI states will also transmit
        X2 <- rgamma(1000000, shape = 2, rate = 2 / aip)
        idx1 <- which(N > 2*X1)
        idx2 <- which(N > 2*(X1+X2))
        beta.m <- sum(idx2) / sum(idx1)
        ## Use this in the calculation of the infection probability
        p.lambda[t, ] <- calc.p.lambda.row(p.beta[t, , ], I.1[t, ] + I.2[t, ], UI.1[t, ] + UI.2[t, ], delta.t, beta.m)
        ## Effective infectious population
        infecting <- I.1[t - 1, ] + I.2[t - 1, ] + (beta.m * (UI.1[t - 1, ] + UI.2[t - 1, ]))
        newly.infecting <- new.I2 + (beta.m * new.UI2)
        ## How many of these can have their contacts followed up?
        indexes <- min(25000 * delta.t * region.population / country.population,
                       theta * sum(newly.infecting))
        ## What is the distribution over the ages
        index.vec <- newly.infecting / sum(newly.infecting) ## multiplied by the indexes?
        ## Use this to get a mean number of secondary contacts
        expect.contacts <- sweep(mixing.model$base.contact.matrices[[which(sapply(mixing.model$intervals, function(x) any(x == t)))]],
                                 2,
                                 S[t - 1, ] + R[t - 1, ],
                                 `*`)
        expect.contacts <- aip * delta.t * sweep(expect.contacts, 1, index.vec, `*`)
        ## Use simulation to see how this mean number of secodary contacts
        ## is impacted by the threshold of only 30 (or 15 contacts to be traced).
        contacts <- rpois(1000000, expect.contacts)
        expect.contacts.trace <- mean(contacts[contacts < chasing.threshold])
        ## Qr[t, ] - expected number of contacts * proportion in R * (1 - prob gaining infection)
        Qr[t, ] <- indexes * expect.contacts.trace * (R[t - 1, ] / (S[t - 1, ] + R[t - 1, ])) * 1
        ## Qs[t, ] - as above
        Qs[t, ] <- indexes * expect.contacts.trace * (S[t - 1, ] / (S[t - 1, ] + R[t - 1, ])) * (1 - p.lambda[t - 1, ])
        ## q[t, ] - probability of being traced given infection at time t
        q[t, ] <- indexes * expect.contacts.trace * (S[t - 1, ] / (S[t - 1, ] + R[t - 1, ])) * p.lambda[t - 1, ]
    }

  ## AGGREGATE THE NUMBER OF NEW INFECTIONS OVER DAYS
                                        #if(!result.by.half.day){
    days <- gl(num.days, intervals.per.day)
    
    if(ncol(NNI) > 1){
        NNI.days <- by(NNI, days, function(x) apply(x, 2, sum))
        NNI.days <- do.call("rbind", NNI.days)
    } else {
        NNI.days <- by(NNI, days, sum)
        NNI.days <- matrix(NNI.days, length(NNI.days), 1)
    }
    
    detach(phi)
    ## }
    seropos <- apply(S, 1, function(x) 1 - (x / population.size))
    seropos <- t(array(seropos, dim = dim(S)[2:1]))[seq(2, nrow(S), by = 2), ]
    list(per.timestep = NNI, per.day = NNI.days, Rt = Rt)
}

SEEIIR.mass.action.model.function <- function(egr, aip, I0.tot, R0.amplitude.param = 0, SEIR.param.list, mixing.model, population.size, intervals.per.day = 2, num.days = 245, start.day = 121, result.by.half.day = F){

  num.age.groups <- length(population.size)
  num.intervals <- num.days * intervals.per.day
  times <- 1:num.intervals
  delta.t <- 1 / intervals.per.day
  
  R.init <- egr * aip * (((egr * SEIR.param.list$latent.period / 2) + 1)^2) / (1 - (1 / (((egr * aip / 2) + 1)^2)))

  R0.amplitude <- R0.amplitude.param * R.init / (1 + cos(2 * pi * (start.day - SEIR.param.list$R0.peak.day) / 365.25))
  
  if(R0.amplitude != 0){
    R0 <- R.init + (R0.amplitude * (cos( 2 * pi * ((times * delta.t) + start.day - SEIR.param.list$R0.peak.day) / 365.25) - cos(2 * pi * (start.day - SEIR.param.list$R0.peak.day) / 365.25)))
  } else R0 <- rep(R.init, length(times))

  ### INITIALISE MODEL
  rho <-  2 * delta.t/ aip
  sigma <- delta.t * 2 / SEIR.param.list$latent.period
  alpha <- exp(egr * delta.t) - 1
  I0 <- I0.tot * mixing.model$initial.age.distribution
  

  I1 <- I0 / (1 + (rho / (alpha + rho)))
  I2 <- I1 * (rho / (alpha + rho))
  E2 <- I1 * ((alpha + rho) / sigma)
  E1 <- E2 * ((alpha + sigma) / sigma)  
  
  NNI <- p.lambda <- S <- I.1 <- I.2 <- E.1 <- E.2 <-R <- matrix(0, num.intervals, num.age.groups)
  
  I.1[1, ] <- I1
  I.2[1, ] <- I2
  E.1[1, ] <- E1
  E.2[1, ] <- E2
  
  S[1, ] <- population.size * SEIR.param.list$susceptibility
  R[1, ] <- population.size * (1 - SEIR.param.list$susceptibility)
  
  p.beta <- array(0, c(num.intervals, num.age.groups, num.age.groups))

  m.max <- min((1:length(mixing.model$intervals))[sapply(1:length(mixing.model$intervals), function(y) num.intervals %in% mixing.model$intervals[[y]])])

  for(m in 1:m.max){
   t.max <- min(num.intervals, mixing.model$intervals[[m]][length(mixing.model$intervals[[m]])])
   t.min <- mixing.model$intervals[[m]][1]
   for(t in t.min:t.max)
   p.beta[t, , ] <- R0[t] * mixing.model$scaled.matrices[[m]] / aip
  }

  p.lambda[1, ] <- calc.p.lambda.row4(p.beta[1, , ], I.1[1, ], I.2[1, ], delta.t, SEIR.param.list$relative.susceptibility.I1.to.I2)
    NNI[1, ] <- p.lambda[1, ] * S[1, ]
  
  for(t in 2:nrow(S)){
    S[t, ] <- S[t - 1, ] * (1 - p.lambda[t -1, ])
    E.1[t, ] <- ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * E.1[t - 1, ]) + (p.lambda[t - 1, ] * S[t - 1, ])
    E.2[t, ] <- ((1 - ((2 * delta.t) / SEIR.param.list$latent.period)) * E.2[t - 1, ]) + ((2 * delta.t / SEIR.param.list$latent.period) * E.1[t - 1, ])
    I.1[t, ] <- (I.1[t - 1, ] * (1 - (2 * delta.t / aip))) + ((2 * delta.t / SEIR.param.list$latent.period) * E.2[t - 1, ])
    I.2[t, ] <- (I.2[t - 1, ] * (1 - (2 * delta.t / aip))) + (I.1[t - 1, ] * 2 * delta.t / aip)
    R[t, ] <- R[t - 1, ] + (I.2[t - 1, ] * 2 * delta.t / aip)
    NNI[t, ] <- p.lambda[t - 1, ] * S[t - 1, ]
    p.lambda[t, ] <- calc.p.lambda.row4(p.beta[t, , ], I.1[t, ], I.2[t, ], delta.t, SEIR.param.list$relative.susceptibility.I1.to.I2)
  }

  ## AGGREGATE THE NUMBER OF NEW INFECTIONS OVER DAYS
  #if(!result.by.half.day){
    days <- gl(num.days, intervals.per.day)
    
  if(ncol(NNI) > 1){
    NNI.days <- by(NNI, days, function(x) apply(x, 2, sum))
    NNI.days <- do.call("rbind", NNI.days)
  } else {
    NNI.days <- by(NNI, days, sum)
    NNI.days <- matrix(NNI.days, length(NNI.days), 1)
  }
  
  seropos <- apply(S, 1, function(x) 1 - (x / population.size))
  seropos <- t(array(seropos, dim = dim(S)[2:1]))[seq(2, nrow(S), by = 2), ]
  list(per.timestep = NNI, per.day = NNI.days, seropos = seropos)

}
