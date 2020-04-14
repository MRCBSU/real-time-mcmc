library(lubridate)

load("mcmc(1).RData")

## Have a look at the posterior range for R0
quantile(posterior.R0, probs = c(0.025, 0.5, 0.975))

## Get Re, the basic reproduction number in the two epochs (pre- and post-lockdown)
## Assuming a fully susceptible population
tlock <- 36
ndays <- dim(NNI[[1]])[2]
posterior.Re <- sapply(1:1000,
                       function(x){
                         c(rep(posterior.R0[10 * x, 1], tlock),
                           rep(params$contact_parameters[10 * x, 2] * posterior.R0[10 * x, 1], ndays - tlock))
                       }
)

## Cumulative incidence as a fraction of the total population
NNI.cum <- apply(NNI[[1]], 3, cumsum) / regions.total.population

## Following two objects should have the same dimension
dim(NNI.cum)
dim(posterior.Re)

## Get a chain of R_t values.
posterior.Rt <- posterior.Re * (1 - NNI.cum)

## Quantiles
qRt <- apply(posterior.Rt, 1, quantile, probs = c(0.025, 0.5, 0.975))

## Get the value of Rt today
dimnames(qRt)[[2]] <- as.character(dates.used)
## Output current R_t
qRt[, which(dimnames(qRt)[[2]] == today())]