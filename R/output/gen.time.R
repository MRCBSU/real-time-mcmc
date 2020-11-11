init.func <- function(vecS, matM){
    if(length(dim(matM)) == 3 && dim(matM)[3] == 1) matM <- as.matrix(matM[,,1])
    if(length(vecS) != nrow(matM)) stop("Dimension mismatch between vecS and matM")
    M.star <- sweep(matM, 1, vecS, `*`)
    eM <- eigen(M.star)$vectors[, 1]
    eM / sum(eM)
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
gen.time.dist <- function(lp, ip, init.vec, M, tol = 0.01, delta.t = 0.5){
    t <- delta.t
    E.1 <- init.vec
    nA <- length(E.1)
    E.2 <- rep(0, nA)
    I.1 <- rep(0, nA)
    I.2 <- rep(0, nA)
    p.lambda <- rep(0, nA)
    while((t * sum(E.1 + E.2 + I.1 + I.2)) > tol){
        E.1.temp <- (1 - ((2 * delta.t) / lp)) * E.1
        E.2.temp <- ((1 - ((2 * delta.t) / lp)) * E.2) + ((2 * delta.t / lp) * E.1)
        I.1.temp <- ((1 - ((2 * delta.t) / ip)) * I.1) + ((2 * delta.t / lp) * E.2)
        I.2.temp <- ((1 - ((2 * delta.t) / ip)) * I.2) + ((2 * delta.t / ip) * I.1)
        E.1 <- E.1.temp
        E.2 <- E.2.temp
        I.1 <- I.1.temp
        I.2 <- I.2.temp
        p.lambda <- rbind(p.lambda,
                          calc.p.lambda.row(M / ip, I.1.temp, I.2.temp, delta.t)
                          )
        t <- t + delta.t
    }
    ## Sum over all the infections (as a proportion of the population size)
    ## gen.time should now be an un-normalised distribution 
    gen.time <- apply(p.lambda, 1, sum)
    gen.time / sum(gen.time)
}
              
