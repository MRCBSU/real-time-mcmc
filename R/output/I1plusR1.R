load("mcmc.RData")

r1.prior.fn <- function(x){
    out <- rep(0, length(x))
    idx <- x>1
    out[idx] <- dgamma(x[idx]-1,shape = 5.5, rate = 1)
    out
}

ld <- density(params$r1_period)
ld$x <- ld$x + 1
pdf("temp.pdf")
plot(r1.prior.fn, xlim = c(0, 20), type = "l", lty = 3, ylim = c(0, max(ld$y)), main = "Prior-posterior for time in R[1]", xlab = "Days")
lines(density(params$r1_period))
dev.off()

testpos <- params$r1_period + params$infectious_period + 3
pdf("temp2.pdf")
plot(density(testpos), main = "Mean total time spent as a prevalent infection", xlab = "Days")
abline(v = mean(testpos), lty =3 , col = "blue")
abline(v = median(testpos), lty = 3, col = "darkgreen")
dev.off()
