suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))

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
  param.to.use <- param[parameter.to.outputs,]
  func.with.param <- function(x) {
    func(x, t(param.to.use))
  }
  dims.to.preserve <- drop.from.names(output.dimnames, c(over, "iteration"))
  mat <- aperm(output, c(over, "iteration", dims.to.preserve))
  return(apply.over.named.array(mat, c(over, "iteration"), func.with.param,
                                target.dimnames))
}

merge.youngest.age.groups <- function(mat, num.to.group = 2, new.name = NULL) {
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
if (!exists("out.dir")) source(file.path(proj.dir, "config.R"))
load(file.path(out.dir, "mcmc.RData"))
rm(dth.dat)
if (!exists("conv")) {
  source(file.path(proj.dir, "R", "output", "gamma_fns.R"))
  source(file.path(proj.dir, "R", "output", "convolution.R"))
}
if (!exists("num.iterations")) source(file.path(proj.dir, "set_up_inputs.R"))
if (!exists("ddelay.mean")) source(file.path(proj.dir, "set_up_pars.R"))

int_iter <- 0:(num.iterations - 1)
## parameter.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.params)
parameter.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= burnin]
## outputs.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.outputs)
outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]
parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)
stopifnot(length(parameter.to.outputs) == length(outputs.iterations)) # Needs to be subset

################################################################
print('Calculating Rt')
Rt.func <- function(vecS, matM){
  if(length(vecS) != nrow(matM)) stop("Dimension mismatch between vecS and matM")
  M.star <- sweep(matM, 2, vecS, `*`)
  max(abs(eigen(M.star, only.values = TRUE)$value))
}
Rt <- list()
R.star <- NULL
M.star <- M <- M.mult <- list()
iterations.for.Rt <- parameter.to.outputs[seq(from = 1, to = length(parameter.to.outputs), length.out = 500)]
## Get the right iterations of the marginal contact parameter chain
m <- params$contact_parameters[iterations.for.Rt, ]
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
  beta.per.region <- ncol(beta) / r
  R0 <- posterior.R0[iterations.for.Rt, , drop = F]
  colnames(R0) <- regions
  for(idir in 1:length(cm.bases)){
	# Read the matrices containing the transmission between age groups (from Ediwn)
    M[[idir]] <- as.matrix(read_tsv(cm.bases[idir], col_names = FALSE, col_types = cols()))
    # Read matrices which select which elements of m to use, add one because R uses 1-based indexing
    # but the model uses 0-based
    M.mult[[idir]] <- as.matrix(read_tsv(cm.mults[idir], col_names = FALSE, col_types = cols())) + 1
  }
  # m.levels[t] is the number of breakpoints passed on day t
  m.levels <- cut(1:ndays, c(0, cm.breaks, Inf))
  beta.levels <- cut(1:ndays, c(0, cm.breaks[-1], Inf))
  names(M) <- names(M.mult) <- NULL
  pop.total <- all.pop[1, ];names(pop.total) <- regions
  for(reg in regions){
    ireg <- which(regions %in% reg)
    for(idir in 1:length(cm.bases)){
	  # M.star[i] <- m[i] * M[i] but select the correct m for the region
      M.star[[idir]] <- array(apply(m, 1,
			function(mm) {
				M[[idir]] * mm[(ireg-1)*m.per.region+M.mult[[idir]]]
			}),
                              dim = c(nA, nA, nrow(m)))
      M.star[[idir]] <- array(apply(beta, 1,
			function(b) {
				M.star[[idir]] * b[(ireg-1)*beta.per.region + max(1, idir - 1)]
			}),
                              dim = c(nA, nA, nrow(m)))
    }
    M.temp <- matrix(M.star[[1]][, , 1, drop = F], dim(M.star[[1]])[1], dim(M.star[[1]])[2])
	# Calculate the scaling R*, the value of Rt for the first matrix (unscaled) TODO: Check this description is correct.
    R.star[ireg] <- Rt.func(regions.total.population[ireg, ] / pop.total[ireg], M.temp)
	# Calculate number of susceptiables for each day
    S <- apply(NNI[[reg]], c(1, 3), cumsum)[,,seq(from = 1, to = length(parameter.to.outputs), length.out = 500)]  ## TxAxI array
    S <- -sweep(S, 2, regions.total.population[ireg, ], `-`) ## TxAxI
	# Calculate the relative Rt values as a function of the next generation matrix for each day
    R.prime <- sapply(1:ndays,
                      function(x) sapply(1:length(iterations.for.Rt), function(i){
                        M.temp <- matrix(M.star[[m.levels[x]]][, , i, drop = F], nA, nA)
                        Rt.func(S[x, , i] / pop.total[ireg], M.temp)
                      }
                      )
    ) ## IxT array
	# Scale the relative Rt values to give the correct R0s
    Rt[[reg]] <- R0[, ireg] * R.prime / R.star[ireg] ## I*T array
  }
  names(R.star) <- regions
}
################################################################
print('Calculating the intrinsic generation time distribution')
colnames(regions.total.population) <- age.labs
rownames(regions.total.population) <- regions
source(file.path(Rfile.loc, "gen.time.R"))
cM <- matrix(M.star[[1]][, , 1, drop = FALSE], nA, nA)
ni <- nrow(R0)
delta.t <- 0.5 ## time-step length
Egt <- Vargt <- list()
for(reg in regions){
    gt <- lapply(1:ni, function(x) gen.time.dist(posterior.lp,
                                                 posterior.aip[iterations.for.Rt[x], ],
                                                 init.func(regions.total.population[reg, ] / pop.total[reg], cM),
                                                 R0[x, reg] * cM / (pop.total[reg] * R.star[reg]))
                 )
    Egt[[reg]] <- sapply(gt, function(gx) sum(gx * seq(delta.t, by = delta.t, length.out = length(gx))))
    Vargt[[reg]] <- sapply(gt, function(gx) sum(gx * (seq(delta.t, by = delta.t, length.out = length(gx))^2)))
    Vargt[[reg]] <- Vargt[[reg]] - (Egt[[reg]]^2)
}
              
################################################################
print('Formatting time series')


## Extract length of dimensions
num.regions <- length(regions)
stopifnot(length(NNI) == num.regions)
num.ages <- dim(NNI[[1]])[1]
stopifnot(num.ages == length(age.labs))
num.days <- dim(NNI[[1]])[2]
num.NNI.iterations <- dim(NNI[[1]])[3]
stopifnot(length(outputs.iterations) == num.NNI.iterations)
output.quantity.dims <- c(num.ages,
                          num.days, num.NNI.iterations, num.regions)
## Calculate all dates used
dates <- seq(
  from = lubridate::as_date(start.date),
  by = 1,
  length = num.days
)
output.dimnames <- list(
  "age" = age.labs,
  "date" = dates,
  "iteration" = outputs.iterations,
  "region" = regions
)

## Get parameters
tbl_params <- as_tibble(params)
tbl_params$iteration <- parameter.iterations

## Get Rt into nice format
Rt.old <- Rt
Rt.dimnames <- output.dimnames[c("date","region")]
Rt.dimnames$iteration <- iterations.for.Rt
Rt.dimnames <- Rt.dimnames[c("iteration", "date","region")]
Rt.dims <- c(dim(Rt.old[[1]]), num.regions)
Rt <- array(
  unlist(Rt.old),
  dim = Rt.dims,
  dimnames = Rt.dimnames
)

## Get NNI into nice format
infections <- array(
  unlist(NNI),
  dim = output.quantity.dims,
  dimnames = output.dimnames
)
rm(NNI)
cum_infections <- infections %>% apply.over.named.array("date", cumsum)

derived.quantity <- function(scaling.param, overdispersion.param, convolution, noise.replicates = 1) {
  ## Calculate deaths
  if (num.ages > 1) {
    raw <- merge.youngest.age.groups(infections)
  } else {
    raw <- infections
  }
  mean <- raw %>%
    apply.convolution(convolution) %>%
    apply.param.to.output(scaling.param, `*`, "age")
  neg.binom.noise <- function(mu, overdispersion, replicates = noise.replicates) {
    stopifnot(is.null(dim(drop(overdispersion))))
    size <- mu / overdispersion
    return(rnbinom(
      n = length(mu) * replicates,
      mu = rep(mu, each = replicates),
      size = rep(size, each = replicates)
    ))
  }
  noise.iterations <- 1:(noise.replicates * length(outputs.iterations))
  noise.dimnames <- dimnames(mean)
  noise.dimnames$iteration <- noise.iterations
  noisy.out <- mean %>%
    apply.param.to.output(overdispersion.param, neg.binom.noise,
                          target.dimnames = noise.dimnames) %>%
    merge.youngest.age.groups(3, "<25")
  return(list(
    mean = mean,
    noisy.out = noisy.out,
    cumulative = noisy.out %>% apply.over.named.array("date", cumsum)
  ))
}
if (hosp.flag == 0) {
  deaths <- noisy_deaths <- cum_deaths <- NULL
} else {
	delay.to.death <- list(
	  incub.mean = 4,
	  incub.sd = 1.41,
	  disease.mean = ddelay.mean,
	  disease.sd = ddelay.sd,
	  report.mean = 0,
	  report.sd = 0
	)
	F.death <- discretised.delay.cdf(delay.to.death, steps.per.day = 1)
  death.data <- derived.quantity(params$prop_case_to_hosp, params$hosp_negbin_overdispersion, F.death)
  deaths <- death.data$mean
  noisy_deaths <- death.data$noisy.out
  cum_deaths <- death.data$cumulative
}
if (gp.flag == 0) {
  hosp <- noisy_hosp <- cum_hosp <- NULL
} else {
	delay.to.hosp <- list(
	  incub.mean = 4,
	  incub.sd = 1.41,
	  disease.mean = hdelay.mean,
	  disease.sd = hdelay.sd,
	  report.mean = 0,
	  report.sd = 0
	)
	F.hosp <- discretised.delay.cdf(delay.to.hosp, steps.per.day = 1)
  hosp.data <- derived.quantity(params$prop_case_to_GP_consultation, params$gp_negbin_overdispersion, F.hosp)
  hosp <- hosp.data$mean
  noisy_hosp <- hosp.data$noisy.out
  cum_hosp <- hosp.data$cumulative
}


## Parse data
print('Loading true data')
load.data <- function(file.names) {
  col.names <- c('date', age.labs)
  names(file.names) <- regions
  to.combine <- dimnames(infections)$age[1:4]
  dat.raw <- suppressMessages(sapply(file.names, read_tsv, col_names = col.names, simplify = FALSE))
  dat.raw[[".id"]] <- "region"
  return(do.call(bind_rows, dat.raw) %>%
    mutate(`<25` = rowSums(.[to.combine])) %>%
    select(-all_of(to.combine)) %>%
    pivot_longer(-c(date, region), names_to = "age")
  )
}
if (hosp.flag == 0) {
  dth.dat <- NULL
} else {
  dth.dat <- load.data(data.files)
}
if (gp.flag == 0) {
  hosp.dat <- NULL
} else {
  hosp.dat <- load.data(gp.data)
}

## Get population
population <- as_tibble(regions.total.population, rownames = "region") %>%
  pivot_longer(-region, names_to = "age")
  
print('Saving results')
save(infections, cum_infections, deaths, cum_deaths, params, dth.dat, noisy_deaths, Rt,
     hosp, noisy_hosp, cum_hosp, population, hosp.dat,
     file = file.path(out.dir, "output_matrices.RData"))
save(Rt, Egt, Vargt, file = file.path(out.dir, "forSPI.RData"))
