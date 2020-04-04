require(coda)
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else if (.Platform$GUI == "RStudio" || Sys.getenv("RSTUDIO") == "1") {
                # We're in RStudio
                return(rstudioapi::getSourceEditorContext()$path)
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

###### WHERE IS THE PROJECT ROUTE DIRECTORY
file.loc <- dirname(thisFile())
proj.dir <- dirname(dirname(file.loc))
source(file.path(proj.dir, "set_up_inputs.R"))

###### WHERE IS THE R FILE DIRECTORY
rfile.dir <- file.loc
source(file.path(rfile.dir, "input_extraction_fns.R"))

###### DIRECTORY CONTAINING MCMC OUTPUT
target.dir <- out.dir

###### HOW IS THE DATA ORGANISED
weekly.data <- FALSE

###### EPIDEMIC DATA TIMING
start.date <- as.Date("17/02/20", format = "%d/%m/%y")
## flusurvey.start.date <- as.Date("25/06/09", format = "%d/%m/%y")

###### Are we accounting for day of the week effects?
day.of.week.effect <- FALSE
if(day.of.week.effect) dow.design.file <- "d_o_w_design_file.txt"
dow.r.breaks <- FALSE

###### Amount of time to be added to posteriors
###### for the average latent and infectious periods
min.waiting.time <- 2

## output type
SMC.output <- FALSE

## run details
a <- nages  ## number of age classes
d <- ndays  ## length in days of the run
i.saved <- 10000 ## number of iterations saved in coda files
i.summary <- 1000 ## number of iterations of summary statistics stored on file

dates.used <- start.date + (0:(d - 1))

r <- length(regions)

regions.total.population <- get.variable.value(target.dir, "regions_population")
## regions.total.population <- sum(regions.total.population) ### TO BE COMMENTED OUT WHEN NOT LOOKING AT ENGLAND ALONE.
## regions.total.population <- regions.total.population[1]

## ## DEFINE PRIOR MODEL ###################################

### WHICH VARIABLES ARE STOCHASTIC?
source(file.path(proj.dir, "set_up_pars.R"))
var.names <- c("exponential_growth_rate_hyper", "l_p_lambda_0_hyper", "prop_susceptible_hyper", "gp_negbin_overdispersion", "hosp_negbin_overdispersion", "latent_period", "infectious_period", "relative_infectiousness", "prop_symptomatic", "contact_parameters", "R0_amplitude_kA", "R0_seasonal_peakday", "exponential_growth_rate", "log_p_lambda_0", "prop_susceptible", "prop_HI_32_to_HI_8", "prop_case_to_GP_consultation", "prop_case_to_hosp", "prop_case_to_death", "importation_rates", "background_GP", "test_sensitivity", "test_specificity", "day_of_week_effects")
### PRIOR INFORMATION
var.priors <- list(distribution = list(NULL, NULL, NULL, list(dgamma), list(dgamma), NULL, list(dgamma), NULL, NULL, list(NULL, dbeta), NULL, NULL, list(dgamma), list(dnorm, dnorm), NULL, NULL, list(dbeta),
                                       list(dbeta), NULL, NULL, NULL, NULL, NULL, NULL), ## informative prior specification
                   parameters = list(NA, NA, NA, pars.eta, pars.eta.h, NA, pars.dI, NA, NA, contact.pars, NA, NA, pars.egr, pars.nu, NA, NA, pars.pgp,
                                     pars.ifr, NA, NA, NA, NA, NA, NA)
                   )
## save the prior specification for use elsewhere.
save(var.names, var.priors, file = file.path(target.dir, "prior.spec.RData"))
## ## ######################################################

## ###### READ IN THE MCMC CHAIN from binary output files
source(file.path(rfile.dir, "readingbinaryfiles.R"))

var.priors <- lapply(var.priors, function(x) x[stochastic.flags])

## p.gp.design.file <- "p_GP_design.txt"
p.gp.design.file <- NULL
p.gp.patch <- FALSE
p.gp.t.breaks <- FALSE
p.gp.r.breaks <- FALSE
## time.NPFS <- c(83)

### FUNCTION FOR ESTABLISHING NUMBER OF PARAMETERS
prior.density <- function(x, distrib, params){
  outval <- 0
  if(identical(deparse(distrib), deparse(dbeta))) outval <- distrib(x, shape1 = params[1], shape2 = params[2])
  if(identical(deparse(distrib), deparse(dgamma))) outval <- distrib(x, shape = params[1], rate = params[2])
  if(identical(deparse(distrib), deparse(dnorm))) outval <- distrib(x, mean = params[1], sd = params[2])
  outval
}

### FUNCTION TO GIVE THE CORRECT NUMBER OF PARAMETERS INVOLVED IN A PRIOR DENSITY
increment.parameters <- function(distrib)
  {
    outval <- 0
    if(identical(deparse(distrib), deparse(dbeta))) outval <- 2
    if(identical(deparse(distrib), deparse(dgamma))) outval <- 2
    if(identical(deparse(distrib), deparse(dnorm))) outval <- 2
    outval
  }

### Possibly useful function
fninvit <- function(p) exp(p) / (1 + exp(p))

### ESTABLISH VARIABLES AS MCMC OBJECTS
for(var.string in var.names[stochastic.flags])
    params[[var.string]] <- as.mcmc(t(params[[var.string]]))

## ## DRAW CODA PLOTS
pdf(file.path(target.dir, "codas.pdf"))
par(mfrow = c(1, 2))
for(inti in 1:npar)
  {
    start.index <- 1
    ## temp.params <- var.priors$parameters[inti]
    for(intj in 1:ncol(params[[inti]]))
      {
        temp.dist <- var.priors$distribution[[inti]][[intj]]
        if(!is.null(temp.dist))
          {
            plot(params[[inti]][, intj], main = paste(parameter.names[inti], intj, sep = "\n"))
            
            ## superimpose the prior densities
            end.index <- start.index + increment.parameters(temp.dist)
            curve(prior.density(x, temp.dist, var.priors$parameters[[inti]][start.index:(end.index - 1)]),
                  min(params[[inti]][, intj]), max(params[[inti]][, intj]), lty = 4, lwd = 1.5, add = TRUE, col = "red"
                  )
          }
        
        start.index <- end.index

      }

  }

## ADD A PLOT FOR THE CHAIN OF R0
R0.func <- function(egr, aip, lp)
  {
    
    out.denom <- 1 - ((((egr * aip) / 2) + 1)^(-2))
    out.numer <- (((egr * lp) / 2) + 1)^2
    
    outval <- egr * aip * out.numer / out.denom
  }

posterior.egr <- params$exponential_growth_rate ## One column for each region, presumably

if(is.null(params$infectious_period))
   {
     posterior.aip <- 2
     prior.aip <- 2
   } else {
     posterior.aip <- min.waiting.time + params$infectious_period   ### Also, one column per region? The average infectious period might have a temporal breakpoint
     prior.aip <- min.waiting.time + rgamma(i.saved, var.priors$parameters$infectious_period[1], rate = var.priors$parameters$infectious_period[2]) ## will need editing for multiple regions      
   }

posterior.lp <- min.waiting.time + value.dl

if(ncol(posterior.aip) != ncol(posterior.egr))
  {
    posterior.R0 <- apply(posterior.egr, 2, R0.func, aip = posterior.aip, lp = posterior.lp)
  } else posterior.R0 <- R0.func(posterior.egr, posterior.aip, posterior.lp)
while(ncol(posterior.R0) < r) posterior.R0 <- cbind(posterior.R0, posterior.R0[, ncol(posterior.R0)])
prior.egr <- rgamma(i.saved, var.priors$parameters$exponential_growth_rate[1], rate = var.priors$parameters$exponential_growth_rate[2])
prior.lp <- min.waiting.time + value.dl
prior.R0 <- R0.func(prior.egr, prior.aip, prior.lp)

for(inti in 1:ncol(posterior.R0))
 {
   main.exp <- parse(text = paste("paste(Trace)~italic(R[0])~paste(,)~paste(", regions[inti], ")", sep = ""))
   plot(posterior.R0[, inti], main = main.exp, type = "l")
   lines(density(prior.R0)$x, density(prior.R0)$y, lty = 4, lwd = 1.5, col = "red")
 }

## ADD A PLOT FOR THE CHAIN OF I0
I0.func <- function(aip, nu, R0, pGP, popn)
    aip * exp(nu) * sum(popn) / (R0 * pGP)
## nu
prior.nu <- rnorm(i.saved, var.priors$parameters$log_p_lambda_0[1], sd = var.priors$parameters$log_p_lambda_0[2])
posterior.nu <- params$log_p_lambda_0 ## Again, presuming one column for each region
## pGP - only the initial propensity in children
if(!is.null(params$prop_case_to_GP_consultation)){
    prior.pGP <- exp(fninvit(rnorm(i.saved, var.priors$parameters$prop_case_to_GP_consultation[1], sd = var.priors$parameters$prop_case_to_GP_consultation[2])))
    nc <- ncol(params$prop_case_to_GP_consultation) 
    posterior.pGP <- exp(fninvit(params$prop_case_to_GP_consultation[, seq(1, nc, by = nc / r)]))
} else {
    prior.pGP <- array(0.1, dim = dim(posterior.R0))
    posterior.pGP <- array(0.1, dim = dim(posterior.R0))
}    ## N

## I0
if(ncol(posterior.aip) == ncol(posterior.nu)){
  posterior.I0 <- t(I0.func(t(posterior.aip), t(posterior.nu), t(posterior.R0), t(posterior.pGP), regions.total.population))
} else {
  posterior.I0 <- array(0, dim = c(nrow(posterior.nu), min(r, ncol(posterior.nu))))
  for(inti in 1:ncol(posterior.nu))
      if(inti <= r)
          posterior.I0[, inti] <- I0.func(posterior.aip, posterior.nu[, inti], posterior.R0[, inti], posterior.pGP[, inti], regions.total.population[inti])
}
prior.I0 <- I0.func(prior.aip, prior.nu, prior.R0, prior.pGP, regions.total.population)
q.prior.I0 <- quantile(prior.I0, probs = c(0.05, 0.95))
prior.I0 <- prior.I0[order(prior.I0)]

for(inti in 1:ncol(posterior.I0))
  {
    main.exp <- parse(text = paste("paste(Trace)~italic(I[0])~paste(,)~paste(", regions[inti], ")", sep = ""))
    plot(posterior.I0[, inti], main = main.exp, type = "l")
    lines(density(prior.I0[floor(i.saved * 0.05):ceiling(i.saved * 0.95)]), lty = 4, lwd = 1.5, col = "red")
  }    

par(mfrow = c(1, 1))

## ADD A PLOT FOR THE LIKELIHOOD CHAIN, SHOULD IT EXIST
if(file.exists(paste(target.dir, "coda_lfx", sep = "")))
  {
    lfx <- as.mcmc(lfx)
    plot(lfx, main = "lfx chain")
  }

dev.off()

## ## ##
## source(paste(rfile.dir, "pGP_curve.R", sep = ""))

## ## ## IF THE PROPENSITY TO CONSULT IS MODELLED
## if(!is.null(p.gp.design.file))
##   {
##     ## DRAW THE PROPENSITY TO CONSULT AS A FUNCTION OF TIME
##     X.mat <- scan(paste(target.dir, p.gp.design.file, sep = ""), )
##     ncs <- ncol(params$prop_case_to_GP_consultation)
##     X.mat <- t(matrix(X.mat, ncs, length(X.mat) / ncs))

##     if(day.of.week.effect){
##       dow.mat <- scan(paste(target.dir, dow.design.file, sep = ""), )
##       ncs <- ncol(params$day_of_week_effects)
##       dow.mat <- t(matrix(dow.mat, ncs, length(dow.mat) / ncs))
##     }
    
##     p.GP <- pGPcurve(params$prop_case_to_GP_consultation,
##                      X.mat,
##                      ifelse(p.gp.r.breaks, r, 1),
##                      time.NPFS,
##                      p.gp.patch,
##                      p.gp.t.breaks,
##                      p.gp.r.breaks,
##                      day.of.week.effect,
##                      params$day_of_week_effects,
##                      dow.mat,
##                      dow.r.breaks,
##                      d)                     
                     
##     xpts <- dates.used[c(1:d, d)] + c(rep(0, d), 1)
    
##     q.p.GP.child <- apply(p.GP$child, 1:2, quantile, probs = c(0.025, 0.5, 0.975))
##     q.p.GP.adult <- apply(p.GP$adult, 1:2, quantile, probs = c(0.025, 0.5, 0.975))

##     for(intr in 1:(dim(q.p.GP.child)[2]))
##       {
##         pdf(paste(target.dir, "p_GP_plots_region", intr, ".pdf", sep = ""))
##         par(mfrow = c(2, 1))
##         plot(xpts, c(q.p.GP.child[3, intr, ], q.p.GP.child[3, intr, dim(q.p.GP.child)[3]]), type = "s", lty = 2, ylab = "P_GP", xlab = "Date", ylim = c(0, 0.55), main = "p_GP Children", lwd = 2)
##         lines(xpts, c(q.p.GP.child[1, intr, ], q.p.GP.child[1, intr, dim(q.p.GP.child)[3]]), type = "s", lwd = 2, lty = 2)
##         lines(xpts, c(q.p.GP.child[2, intr, ], q.p.GP.child[2, intr, dim(q.p.GP.child)[3]]), type = "s", lwd = 2.5)
##         plot(xpts, c(q.p.GP.adult[3, intr, ], q.p.GP.adult[3, intr, dim(q.p.GP.adult)[3]]), type = "s", lty = 2, ylab = "P_GP", xlab = "Date", ylim = c(0, 0.55), main = "p_GP Adults", lwd = 2)
##         lines(xpts, c(q.p.GP.adult[1, intr, ], q.p.GP.adult[1, intr, dim(q.p.GP.adult)[3]]), type = "s", lty = 2, lwd = 2)
##         lines(xpts, c(q.p.GP.adult[2, intr, ], q.p.GP.adult[2, intr, ncol(q.p.GP.adult)]), type = "s", lwd = 2.5)
        
##         load(paste(proj.dir, "Data/FluSurvey/prop_GP.RData", sep = ""))
      
##         flusurvey.points <- flusurvey.start.date + (7 * seq(0:26))
        
##         points(flusurvey.points, pGP[1:length(flusurvey.points)], col = "red", pch = "+")
##         lines(flusurvey.points, pGP[1:length(flusurvey.points)], col = "red")
##         points(flusurvey.points, vGP[1:length(flusurvey.points)], col = "blue", pch = "+")
##         lines(flusurvey.points, vGP[1:length(flusurvey.points)], col = "blue")
        
        
##         dev.off()
##       }
##   } else { ### fix-ups to make the other bits of R code work when there's no GP data
    
##     p.GP.child <- matrix(0.1, 1, 10000) ### important for calculating I0 check mod_pars.txt
    
##   }

save.image(file.path(target.dir, "mcmc.RData"))
