gamma.entropy <- function(a)    a[1] - log(a[2]) + lgamma(a[1]) + ((1 - a[1]) * digamma(a[1]))

gamma.mean.var.params <- function(a) c(a[1] / a[2], a[1] / (a[2]^2))

gamma.shape.rate.params <- function(mean.var) c((mean.var[1]^2) / mean.var[2], mean.var[1] / mean.var[2])

gamma.approx.kb <- function(x, smallf){
i <- 0:(length(smallf) - 1)
approx.smallf <- sapply(i, function(y, gamma.params) integrate(dgamma, lower = max(y - 0.5, 0), upper = ifelse(y == max(i), Inf, y + 0.5), shape = gamma.params[1], rate = gamma.params[2])$value, gamma.params = x)
kb.div <- smallf * log(smallf / approx.smallf)
sum(kb.div)
}
