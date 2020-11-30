eval.drop <- function(popn, mat.before, mat.after){
    ngm.before <- popn * mat.before
    ngm.after <- popn * mat.after
    eigen(ngm.after, only.values = TRUE)$values[1] / eigen(ngm.before, only.values = TRUE)$values[1]
}

R0.drop <- lapply(mat.intervention,
                  function(xmat) apply(susceptible_popns,
                                       which(names(dim(susceptible_popns)) %in% c("iteration", "region")),
                                       eval.drop, mat.before, xmat)
                  )

lapply(R0.drop, function(x) apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975)))
