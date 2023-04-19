### File originally formed part of tidy_output.R but sections could be made re-usable for purposes of scenario generation
Rt.func <- function(vecS, matM){
  if(length(vecS) != nrow(matM)) stop("Dimension mismatch between vecS and matM")
  M.star <- sweep(matM, 1, vecS, `*`)
  max(abs(eigen(M.star, only.values = TRUE)$value))
}
colnames(regions.total.population) <- age.labs
rownames(regions.total.population) <- regions
Rt <- list()
R.star <- vector("list", nr)
gt <- vector("list", nr)
M.star <- M <- M.mult <- list()
iterations.for.Rt <- parameter.to.outputs[seq(from = 1, to = length(parameter.to.outputs), length.out = 500)]
outputs.for.Rt <- which(parameter.to.outputs %in% iterations.for.Rt)
## Get the right iterations of the marginal contact parameter chain
m <- params$contact_parameters[iterations.for.Rt, ]
if(beta.update)
    beta <- exp(params$log_beta_rw[iterations.for.Rt, ] %*% t(beta.design))
## Multiply by the design matrix if applicable
if(rw.flag)
    m <- m %*% t(m.design)
## Inverse transformation
m <- exp(m)
# m now has the parameters on natural scale
if(ncol(m) %% r != 0) {
  warning('Number of m parameters is not a multiple of number of regions, cannot caclulate Rt')
} else {
  m.per.region <- ncol(m) / r
  if(beta.update) beta.per.region <- ncol(beta) / r
  R0 <- posterior.R0[iterations.for.Rt, , drop = F]
  ni <- nrow(R0)
  colnames(R0) <- regions
  for(idir in 1:length(cm.bases)){
	# Read the matrices containing the transmission between age groups (from Edwin)
    M[[idir]] <- as.matrix(read_tsv(cm.bases[idir], col_names = FALSE, col_types = cols()))
    # Read matrices which select which elements of m to use, add one because R uses 1-based indexing
    # but the model uses 0-based
    M.mult[[idir]] <- as.matrix(read_tsv(cm.mults[idir], col_names = FALSE, col_types = cols())) + 1
  }
  # m.levels[t] is the number of breakpoints passed on day t
  m.levels <- cut(1:ndays, c(0, cm.breaks, Inf))
  if(beta.update) {
      beta.levels <- cut(1:ndays, c(0, beta.breaks, Inf))
      beta.map <- as.numeric(cut(1:length(cm.bases),c(0, which(cm.breaks %in% beta.breaks), Inf)))
  }
  names(M) <- names(M.mult) <- NULL
  pop.total <- all.pop[1, ];names(pop.total) <- regions
  for(reg in regions){
    ireg <- which(regions %in% reg)
    for(idir in 1:length(cm.bases)){
        ## M.star[i] <- m[i] * M[i] but select the correct m for the region
        M.star[[idir]] <- array(apply(m, 1,
                                      function(mm) {
                                          M[[idir]] * mm[(ireg-1)*m.per.region+M.mult[[idir]]]
                                      }),
                                dim = c(nA, nA, nrow(m)))
        if(beta.update)
            M.star[[idir]] <- array(sapply(1:nrow(beta),
                                           function(b) {
                                               M.star[[idir]][,,b] * beta[b,(ireg-1)*beta.per.region + beta.map[idir]]
                                           }),
                                    dim = c(nA, nA, nrow(m)))
    }
    ## M.temp <- matrix(M.star[[1]][, , 1, drop = F], dim(M.star[[1]])[1], dim(M.star[[1]])[2])
    ## ## Calculate the scaling R*, the value of Rt for the first matrix (unscaled) TODO: Check this description is correct.
    ## R.star[ireg] <- Rt.func(regions.total.population[ireg, ] / pop.total[ireg], M.temp)
    R.star[[ireg]] <- apply(M.star[[1]], 3, function(x) Rt.func(regions.total.population[ireg, ] / pop.total[ireg], x))
    ## Generation time distribution based on the initial contact matrix
    gt[[ireg]] <- lapply(1:ni, function(x) gen.time.dist(posterior.lp,
                                                         posterior.aip[iterations.for.Rt[x], ],
                                                         init.func(regions.total.population[reg,] / pop.total[reg], M.star[[1]][,,x]),
                                                         R0[x,reg] * M.star[[1]][,,x] / (pop.total[reg] * R.star[[ireg]][x]))
                         )
    ## Calculate number of susceptibles for each day
    ## S <- apply(NNI[[reg]][,,outputs.for.Rt, drop = FALSE], c(1, 3), cumsum)  ## TxAxI array
    ## S <- -sweep(S, 2, regions.total.population[ireg, ], `-`) ## TxAxI
    S <- (regions.total.population[ireg, ] * (1 - sero[[reg]][,,outputs.for.Rt,drop=FALSE])) %>%
        aperm(c(2, 1, 3))
    ## Calculate the relative Rt values as a function of the next generation matrix for each day
    R.prime <- sapply(1:ndays,
                      function(x) sapply(1:length(iterations.for.Rt),
                                         function(i){
                                             MM <- matrix(M.star[[m.levels[x]]][, , i, drop = F], nA, nA)
                                             Rt.func(S[x, , i] / pop.total[ireg], MM)
                                         }
                                         )
                      ) ## IxT array
    ## Scale the relative Rt values to give the correct R0s
    Rt[[reg]] <- R0[, ireg] * R.prime / R.star[[ireg]] ## I*T array
  }
  names(R.star) <- regions
}
