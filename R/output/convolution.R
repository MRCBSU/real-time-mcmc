## Load in convolution functions for estimating Hospitalisations, ICUs and Deaths
dyn.load(file.path(proj.dir, "C-exts", "convolve", "convolve.so"))
## Define a convolution function (coded in C++)
conv <- function(a, b)
    .C("convolute",
       as.double(a),
       as.integer(length(a)),
       as.double(b),
       as.integer(length(b)),
       ab = double(length(a) + length(b) - 1))$ab

### FUNCTION REQUIRED BY reporting.model
discretised.delay.cdf <- function(delay.list, steps.per.day = 2, zero.prob = 1e-5){
  
### get overall delay mean
  overall.delay.mean <- do.call("sum",
                                delay.list[grep("mean$", names(delay.list))]
                                )
  
### get overall delay variance
  overall.delay.var <- do.call(function(...) sum(c(...)^2),
                               delay.list[grep("sd$", names(delay.list))]
                               )
  
### get gamma shape rate parameters
  shape.rate <- gamma.shape.rate.params(c(overall.delay.mean, overall.delay.var))
  
### get bound on how many intervals we realistically need to consider
  effective.upper.limit <- ceiling(qgamma(1 - zero.prob, shape = shape.rate[1], rate = shape.rate[2])) * steps.per.day
  
### get discretised cdf
  pgamma(c((1:effective.upper.limit) / steps.per.day, Inf), shape = shape.rate[1], rate = shape.rate[2]) - pgamma((0:effective.upper.limit) / steps.per.day, shape = shape.rate[1], rate = shape.rate[2])
  
}
