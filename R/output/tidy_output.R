suppressMessages(library(tidyverse))

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
if(!exists("proj.dir")){
  file.loc <- dirname(thisFile())
  proj.dir <- dirname(dirname(file.loc))
}
if (!exists("out.dir")) source(file.path(proj.dir, "config.R"))
load(file.path(out.dir, "mcmc.RData"))
if (!exists("conv")) {
  source(file.path(proj.dir, "R", "output", "gamma_fns.R"))
  source(file.path(proj.dir, "R", "output", "convolution.R"))
}
if (!exists("num.iterations")) source(file.path(proj.dir, "set_up_inputs.R"))
if (!exists("ddelay.mean")) source(file.path(proj.dir, "set_up_pars.R"))


parameter.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.params)
outputs.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.outputs*2)
parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)
stopifnot(length(parameter.to.outputs) == length(outputs.iterations)) # Needs to be subset

################################################################
Rt.func <- function(vecS, matM){
  if(length(vecS) != nrow(matM)) stop("Dimension mismatch between vecS and matM")
  M.star <- sweep(matM, 2, vecS, `*`)
  max(abs(eigen(M.star, only.values = TRUE)$value))
}
Rt <- list()
M.star <- M <- M.mult <- list()
m <- params$contact_parameters[parameter.to.outputs, ]
R0 <- posterior.R0[parameter.to.outputs, , drop = F]
for(idir in 1:length(cm.bases)){
  M[[idir]] <- as.matrix(read_tsv(cm.bases[idir], col_names = FALSE))
  M.mult[[idir]] <- as.matrix(read_tsv(cm.mults[idir], col_names = FALSE)) + 1
  M.star[[idir]] <- array(apply(m, 1, function(mm) M[[idir]] * mm[M.mult[[idir]]]),
                          dim = c(nA, nA, nrow(m)))
}
m.levels <- cut(1:ndays, c(0, cm.breaks, Inf))
names(M) <- names(M.mult) <- names(M.star)
pop.total <- all.pop[1, ]
for(reg in regions){
  ireg <- which(regions %in% reg)
  M.temp <- matrix(M.star[[1]][, , 1, drop = F], dim(M.star[[1]])[1], dim(M.star[[1]])[2])
  R.star <- Rt.func(regions.total.population[ireg, ] / pop.total[ireg], M.temp)
  S <- apply(NNI[[reg]], c(1, 3), cumsum)  ## TxAxI array
  S <- -sweep(S, 2, regions.total.population[ireg, ], `-`) ## TxAxI
  R.prime <- sapply(1:ndays,
                    function(x) sapply(1:i.summary, function(i){
                      M.temp <- matrix(M.star[[m.levels[x]]][, , i, drop = F], nA, nA)
                      Rt.func(S[x, , i] / pop.total[ireg], M.temp)
                    }
                    )
  ) ## IxT array
  Rt[[reg]] <- R0[, ireg] * R.prime / R.star ## I*T array
}
################################################################

## Constants and settings

delay.to.death <- list(
  incub.mean = 4,
  incub.sd = 1.41,
  disease.mean = ddelay.mean,
  disease.sd = ddelay.sd,
  report.mean = rdelay.mean,
  report.sd = rdelay.sd
)
F.death <- discretised.delay.cdf(delay.to.death, steps.per.day = 1)

## Extract length of dimensions
NNI <- NNI[,,1:length(outputs.iterations)]
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
Rt.dimnames <- output.dimnames[c("iteration","date","region")]
Rt.dims <- sapply(output.dimnames[c("iteration","date","region")], length)
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
cum_infections <- infections %>% apply.over.named.array("date", cumsum)

## Calculate deaths
if (num.ages > 1) {
  deaths <- merge.youngest.age.groups(infections)
} else {
  deaths <- infections
}
deaths <- deaths %>%
  apply.convolution(F.death) %>%
  apply.param.to.output(params$prop_case_to_hosp, `*`, "age")
noise.replicates <- 5
neg.binom.noise <- function(mu, overdispersion, replicates = noise.replicates) {
  stopifnot(is.null(dim(drop(overdispersion))))
  size <- mu / params$hosp_negbin_overdispersion
  return(rnbinom(
    n = length(mu) * replicates,
    mu = rep(mu, each = replicates),
    size = rep(size, each = replicates)
  ))
}
noise.iterations <- 1:(noise.replicates * length(outputs.iterations))
noise.dimnames <- dimnames(deaths)
noise.dimnames$iteration <- noise.iterations
noisy_deaths <- deaths %>%
  apply.param.to.output(params$hosp_negbin_overdispersion, neg.binom.noise,
                        target.dimnames = noise.dimnames) %>%
  merge.youngest.age.groups(3, "<25")
cum_deaths <- noisy_deaths %>% apply.over.named.array("date", cumsum)

## Parse data
dth.col.names <- c('date', age.labs)
names(data.files) <- regions
to.combine <- dimnames(infections)$age[1:4]
dth.dat.raw <- suppressMessages(sapply(data.files, read_tsv, col_names = dth.col.names, simplify = FALSE))
dth.dat.raw[[".id"]] <- "region"
dth.dat <- do.call(bind_rows, dth.dat.raw) %>%
  mutate(`<25` = rowSums(.[to.combine])) %>%
  select(-to.combine) %>%
  pivot_longer(-c(date, region), names_to = "age")

## Get population
colnames(regions.total.population) <- age.labs
rownames(regions.total.population) <- regions
population <- as_tibble(regions.total.population, rownames = "region") %>%
  pivot_longer(-region, names_to = "age")
  
save(infections, cum_infections, deaths, cum_deaths, params, dth.dat, noisy_deaths, Rt,
     population,
     file = file.path(out.dir, "output_matrices.RData"))
