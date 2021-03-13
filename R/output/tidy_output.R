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
## save.image("tmptmp.RData")
################################################################
print('Calculating Rt')
source(file.path(Rfile.loc, "gen.time.R"))
source(file.path(Rfile.loc, "calculating_Rt.R"))
################################################################
print('Calculating the intrinsic generation time distribution')
delta.t <- 0.5 ## time-step length
Egt <- Vargt <- list()
names(gt) <- regions
for(reg in regions){
    ## gt <- lapply(1:ni, function(x) gen.time.dist(posterior.lp,
    ##                                              posterior.aip[iterations.for.Rt[x], ],
    ##                                              init.func(regions.total.population[reg, ] / pop.total[reg], M.star[[1]][,,x]),
    ##                                              R0[x, reg] * M.star[[1]][,,x] / (pop.total[reg] * R.star[[reg]][x]))
    ##              )
    Egt[[reg]] <- sapply(gt[[reg]], function(gx) sum(gx * seq(delta.t, by = delta.t, length.out = length(gx))))
    Vargt[[reg]] <- sapply(gt[[reg]], function(gx) sum(gx * (seq(delta.t, by = delta.t, length.out = length(gx))^2)))
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
nice.array <- function(x)
    array(
        unlist(x),
        dim = output.quantity.dims,
        dimnames = output.dimnames)
infections <- nice.array(NNI); rm(NNI)
sero <- nice.array(sero)
if(vacc.flag) vacc.infections <- nice.array(DNNI); rm(DNNI)
cum_infections <- infections %>% apply.over.named.array("date", cumsum)
if(dths.flag){
    deaths2 <- nice.array(Deaths)
    rm(Deaths)
}
if(prev.flag){
    prevalence <- nice.array(Prev)
    rm(Prev)
	if (!exists("prev.dat")) {
	  if (!exists("prev.dat.file")) prev.dat.file <- paste0(data.dirs["prev"], "/", date.prev, "_", last.prev.day, "_ons_dat.csv")
	  prev.dat <- read_csv(prev.dat.file)
	}
} else prev.dat <- NULL
derived.quantity <- function(scaling.param, scaling.idxs = c("date", "age"),
                             overdispersion.param, convolution,
                             series = infections, dow = NULL,
                             noise.replicates = 1, observe.babies = FALSE, merge.youngest.grps = 3, merge.youngest.label = "<25") {
    if(!observe.babies & num.ages > 1){
        ## ## Merging of youngest age groups needs to be done prior to adding negbin noise.
        raw <- merge.youngest.age.groups(series)
    } else {
        raw <- series
    }
    mean <- raw %>%
        apply.convolution(convolution) %>%
        apply.param.to.output(scaling.param, `*`, scaling.idxs)
    if(!is.null(dow))
        mean <- apply.param.to.output(mean, dow, `*`, "date")
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
        apply.param.to.output(array(overdispersion.param, dim = dim(overdispersion.param)),
                              neg.binom.noise,
                              target.dimnames = noise.dimnames) %>%
        merge.youngest.age.groups(merge.youngest.grps, merge.youngest.label)
    return(list(
        mean = mean,
        noisy.out = noisy.out,
        cumulative = noisy.out %>% apply.over.named.array("date", cumsum)
    ))
}

pad.fullrange <- function(in.vec, breaks, rmax){
    lths <- diff(c(0, breaks, rmax))
    unlist(lapply(1:(length(breaks)+1), function(y) rep(in.vec[y], lths[y])))
    }

if (!dths.flag) {
  deaths <- noisy_deaths <- cum_deaths <- NULL
} else {
    delay.to.death <- list(
        incub.mean = 4,
        incub.sd = 1.41,
        disease.mean = ddelay.mean,
        disease.sd = ddelay.sd,
        report.mean = rdelay.mean,
        report.sd = rdelay.sd
    )
    F.death <- discretised.delay.cdf(delay.to.death, steps.per.day = 1)
    ## Do we have a single IFR for each region, or do we need a design matrix?
    if(single.ifr){
        ifr <- array(params$prop_case_to_hosp, dim = c(dim(params$prop_case_to_hosp), ndays))
    } else {
        ## the manipulation of the IFR below was originally written in tidy syntax but it caused memory problems
        ifr <- expit(params$prop_case_to_hosp %*% t(model.matrix(lm.TA2)))
        ## ncols gives the number of different values for the IFR, need to spread them out over time and age
        ifr <- array(ifr, dim = c(dim(ifr)[1], nA - 1, length(tbreaks.ifr) + 1))
        ## expand the date dimension according to the specified breakpoints
        ifr <- apply(ifr, 1:2, pad.fullrange, breaks = tbreaks.ifr, rmax = ndays)
        ifr <- aperm(ifr, c(2, 3, 1))
    }
    merge.flg <- grepl("adjusted", data.desc)
    death.data <- derived.quantity(ifr,
                                   overdispersion.param = params$hosp_negbin_overdispersion,
                                   convolution = F.death,
                                   merge.youngest.grps = 3 + merge.flg,
                                   merge.youngest.label = ifelse(merge.flg, "<45", "<25"))
    deaths <- death.data$mean
    noisy_deaths <- death.data$noisy.out
    cum_deaths <- death.data$cumulative
}
if (!cases.flag) {
  cases <- noisy_case <- cum_cases <- NULL
} else {
    delay.to.case <- list(
        incub.mean = 4,
        incub.sd = 1.41,
        disease.mean = 0,
        disease.sd = 0,
        report.mean = ldelay.mean,
        report.sd = ldelay.sd
    )
    F.case <- discretised.delay.cdf(delay.to.case, steps.per.day = 1)
    if(nA == 1){
        icr <- array(params$prop_case_to_GP_consultation, dim = c(dim(params$prop_case_to_GP_consultation), nA))
    } else {
        icr <- expit(params$prop_case_to_GP_consultation %*% t(model.matrix(ex6)))
        ## icr has regional and age variation at present
        icr <- array(icr, dim = c(dim(icr)[1], dim(icr)[2] / r, r))
        icr <- aperm(icr, c(1, 3, 2))
        if(gp.flag & nA != 1){ ## this will usually mean the are age breakpoints
            icr <- apply(icr, 1:2, pad.fullrange, breaks = abreaks.icr, nA)
            icr <- aperm(icr, c(2, 1, 3))
        }
    }
    dow <- exp(params$day_of_week_effects %*% t(lm.mat))
    dow <- array(dow, dim = c(dim(dow)[1], end.gp - start.gp + 1))
    dow <- apply(dow, 1, pad.fullrange, breaks = start.gp:(end.gp - 1), ndays)
    case.data <- derived.quantity(icr,
                                  scaling.idxs = c("region", "age"),
                                  overdispersion.param = params$gp_negbin_overdispersion,
                                  convolution = F.case,
                                  dow = t(dow),
                                  observe.babies = TRUE,
                                  merge.youngest.label = "<15")
    case <- case.data$mean
    noisy_case <- case.data$noisy.out
    cum_case <- case.data$cumulative
}


## Parse data
print('Loading true data')
load.data <- function(file.names, idx.age.to.combine = 1:4, label.age.to.combine = "<25") {
  col.names <- c('date', age.labs)
  names(file.names) <- regions
  dat.raw <- suppressMessages(sapply(file.names, read_tsv, col_names = col.names, simplify = FALSE))
  dat.raw[[".id"]] <- "region"
  tbl_dat <- do.call(bind_rows, dat.raw)
  if (num.ages > 1) {
    to.combine <- dimnames(infections)$age[idx.age.to.combine]
    tbl_dat <- tbl_dat %>%
      mutate(!!label.age.to.combine := rowSums(.[to.combine])) %>%
      select(-all_of(to.combine))
  }
  return(tbl_dat %>% pivot_longer(-c(date, region), names_to = "age"))
}
if (hosp.flag == 0 || merge.flg) {
  dth.dat <- NULL
} else  {
  dth.dat <- load.data(data.files)
}
if (gp.flag == 0) {
  case.dat <- NULL
} else {
  case.dat <- load.data(cases.files, 1:3, "<15")
}

## Get population
population <- as_tibble(regions.total.population, rownames = "region") %>%
  pivot_longer(-region, names_to = "age")
  
print('Saving results')
save.objs <- c("infections", "cum_infections", "vacc.infections", "deaths", "cum_deaths", "prevalence", "params", "dth.dat", "noisy_deaths", "Rt",
     "case", "noisy_case", "cum_case", "population", "case.dat", "ifr", "prev.dat")
save.exists <- save.objs[sapply(save.objs, exists)]
save(list = save.exists,
     file = file.path(out.dir, "output_matrices.RData"))
save(Rt, Egt, Vargt, file = file.path(out.dir, "forSPI.RData"))

uneeded <- setdiff(save.exists, ls())
rm(uneeded)

