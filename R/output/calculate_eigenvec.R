suppressMessages(library(tidyverse))
suppressPackageStartupMessages(suppressWarnings(require(cubelyr)))
suppressMessages(library(Matrix))
suppressMessages(extract <- R.utils::extract)

## Location of this script
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

###############################################################

## Functions
logit <- function(p) log(p/(1-p))
expit <- function(x) exp(x) / (1 + exp(x))

drop.from.names <- function(x, value) {
  names(x)[!(names(x) %in% value)]
}

apply.over.named.array <- function(arr, over, func,
                                   target.dimnames = dimnames(arr)) {
  dims.to.preserve <- drop.from.names(dimnames(arr), over)
  new.dim.order <- c(over, dims.to.preserve)
  target.dimnames <- target.dimnames[new.dim.order]
  target.dims <- sapply(target.dimnames, length)
  return(
    apply(arr, dims.to.preserve, func) %>%
      array(dim = target.dims, dimnames = target.dimnames)
  )
}

apply.param.to.output <- function(output, param, func, over = NULL,
                                  target.dimnames = dimnames(output)) {
  param.to.use <- extract(param, "1" = parameter.to.outputs, drop = FALSE)
  func.with.param <- function(x) {
    func(x, aperm(param.to.use, perm = length(dim(param.to.use)):1))
  }
  dims.to.preserve <- drop.from.names(output.dimnames, c(over, "iteration"))
  mat <- aperm(output, c(over, "iteration", dims.to.preserve))
  return(apply.over.named.array(mat, c(over, "iteration"), func.with.param,
                                target.dimnames))
}

merge.youngest.age.groups <- function(mat, num.to.group = 2, new.name = NULL) {
  if (num.ages < num.to.group) {return(mat)}
  add.function <- function(x) {
    to.merge <- x[1:num.to.group]
    to.preserve <- x[(num.to.group+1):length(x)]
    return(c(sum(to.merge), to.preserve))
  }
  dims.to.preserve <- drop.from.names(output.dimnames, "age")
  old.age.group.names <- dimnames(mat)$age
  if (is.null(new.name)){
    new.name <- paste(
      old.age.group.names[1:num.to.group], 
      collapse = ","
    )
  }
  result <- apply(mat, dims.to.preserve, add.function)
  dimnames(result)[[1]][1] <- new.name
  names(dimnames(result))[1] <- "age"
  return(result)
}

apply.convolution <- function(start, func, over = "date") {
  result.length <- length(dimnames(start)[[over]])
  name.over <- dimnames(start)[[over]]
  dims.to.preserve <- drop.from.names(output.dimnames, over)
  result <- start %>%
    apply(dims.to.preserve, conv, b = func)
  result <- result[1:result.length,,,,drop = FALSE]
  dimnames(result)[[1]] <- name.over
  names(dimnames(result))[1] <- over
  return(result)
}

###############################################################

## Load files
print('Loading results')
if(!exists("proj.dir")){
  file.loc <- dirname(thisFile())
  proj.dir <- dirname(dirname(file.loc))
}

load(file.path(out.dir, "mcmc.RData"))
rm(dth.dat)
if (!exists("conv")) {
  source(file.path(proj.dir, "R", "output", "gamma_fns.R"))
  source(file.path(proj.dir, "R", "output", "convolution.R"))
}


int_iter <- 0:(num.iterations - 1)
## parameter.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.params)
parameter.iterations.raw <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= burnin]
parameter.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= min.iteration]
## outputs.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.outputs)
outputs.iterations.raw <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]
outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= min.iteration]
parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)
stopifnot(length(parameter.to.outputs) == length(outputs.iterations)) # Needs to be subset
## save.image("tmptmp.RData")
################################################################
print('Calculating Rt')
source(file.path(Rfile.loc, "gen.time.R"))

### File originally formed part of tidy_output.R but sections could be made re-usable for purposes of scenario generation
Rt.func <- function(vecS, matM){
  if(length(vecS) != nrow(matM)) stop("Dimension mismatch between vecS and matM")
  M.star <- sweep(matM, 1, vecS, `*`)
  max(abs(eigen(M.star, only.values = TRUE)$value))
}
eigenvector.func <- function(vecS, matM){
  if(length(vecS) != nrow(matM)) stop("Dimension mismatch between vecS and matM")
  M.star <- sweep(matM, 1, vecS, `*`)
  eigen_mat <- eigen(M.star, only.values = FALSE)
  eigen_mat$vectors[,which.max(abs(eigen_mat$value))]
}
colnames(regions.total.population) <- age.labs
rownames(regions.total.population) <- regions
Rt <- list()
R.eigenvec <- list()
Contact.eigenvec <- list()
Contact.eigenval <- list()
alt_Contact.eigenvec <- list()
alt_Contact.eigenval <- list()
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
    M.mult[[idir]] <- as.matrix(read_tsv(
      # cm.mults[idir]
      "/home/phe.gov.uk/joel.kandiah/mcmc/real-time-mcmc/contact_mats/ag8_mult_mod4levels0.txt", col_names = FALSE, col_types = cols())) + 1
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
    R.eigenvec[[reg]] <- sapply(1:ndays,
                      function(x) sapply(1:length(iterations.for.Rt),
                                         function(i){
                                             MM <- matrix(M.star[[m.levels[x]]][, , i, drop = F], nA, nA)
                                             eigenvector.func(S[x, , i] / pop.total[ireg], MM)
                                         }
                                    )
                      )

    Contact.eigenvec[[reg]] <- sapply(1:ndays,
                                            function(x) {
                                             a <- eigen(matrix(M[[m.levels[x]]],nA, nA))
                                             a$vectors[,which.max(abs(a$values))]
                                         }
                                         )
    
    Contact.eigenval[[reg]] <- sapply(1:ndays,
                                            function(x) {
                                             a <- eigen(matrix(M[[m.levels[x]]],nA, nA), only.values = TRUE)
                                             max(abs(a$values))
                                         }
                                         )
    
    alt_Contact.eigenvec[[reg]] <- sapply(1:ndays,
                                            function(x) {
                                             a <- eigen(sweep(matrix(M[[m.levels[x]]],nA, nA), 1, regions.total.population[ireg, ], `*`))
                                             a$vectors[,which.max(abs(a$values))]
                                         }
                                         )
    
    alt_Contact.eigenval[[reg]] <- sapply(1:ndays,
                                            function(x) {
                                             a <- eigen(sweep(matrix(M[[m.levels[x]]],nA, nA), 1, regions.total.population[ireg, ], `*`), only.values = TRUE)
                                             max(abs(a$values))
                                         }
                                         )

    ## Scale the relative Rt values to give the correct R0s
    Rt[[reg]] <- R0[, ireg] * R.prime / R.star[[ireg]] ## I*T array
  }
#   names(R.star) <- regions
}

save(Rt, R.eigenvec, Contact.eigenvec, Contact.eigenval, alt_Contact.eigenvec, alt_Contact.eigenval, ndays, nA, nr, iterations.for.Rt, file = "eigenvec.RData")
