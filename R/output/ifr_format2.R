## Written for inclusion in report-updated.Rmd
require(tidyverse)

## Presumes tmp.RData mcmc.RData and output_matrices.RData are all loaded in to the workspace.
if(!exists("out.dir"))
    load("tmp.RData") ## There is a function somewhere in the repo that will guess sensible (relative) file paths. Update to use this.
if(!exists("ifr"))
    load(file.path(out.dir, "output_matrices.RData"))

load(file.path(out.dir, "output_matrices.RData"))
load(file.path(out.dir, "mcmc.RData"))

## ## Base R proved quicker than tidyr equivalents, to...
## Aggregate infections over region
## Aggegrate the two youngest age groups
drop.from.names <- function(x, value) {
  names(x)[!(names(x) %in% value)]
}
merge.youngest.age.groups <- function(mat, output.dimnames, num.to.group = 2, new.name = NULL) {
  if (nA < num.to.group) {return(mat)}
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
sum.up.and.merge <- function(inf){
    inf_by_age <- apply(inf, c("age", "date", "iteration"), sum)
    inf_by_age <- merge.youngest.age.groups(inf_by_age, output.dimnames = dimnames(inf_by_age))
    names(dimnames(inf_by_age))[1] <- "age"
    dimnames(inf_by_age)$age <- dimnames(ifr)$age
    inf_by_age
}

if(!exists("parameter.to.outputs")){
    int_iter <- 0:(num.iterations - 1)
    ## parameter.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.params)
    parameter.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.params)) & int_iter >= burnin]
    ## outputs.iterations <- seq(from = burnin, to = num.iterations-1, by = thin.outputs)
    outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]
    parameter.to.outputs <- which(parameter.iterations %in% outputs.iterations)
}
if(is.null(names(dim(ifr)))){
    names(dim(ifr)) <- c("iteration", "age", "date")
    dimnames(ifr) <- list(iteration = parameter.iterations, age = dimnames(deaths)$age, date = dimnames(deaths)$date)
}

outputs.iterations <- int_iter[(!((int_iter + 1 - burnin) %% thin.outputs)) & int_iter >= burnin]
num.days <- dim(NNI[[1]])[2]
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

## Only look at IFR iterations for which the infections array has been evaluated
ifr.tab <- R.utils::extract(ifr, indices=list(parameter.to.outputs), dims="iteration")

## Format infections
inf_by_age <- sum.up.and.merge(infections) %>% aperm(names(dimnames(ifr.tab)))
if(vacc.flag){
	dinf_by_age <- sum.up.and.merge(vacc.infections) %>% aperm(names(dimnames(ifr.tab)))
	## Calculate a scaled version of the IFR
	ifr.tab <- ifr.tab * dinf_by_age / inf_by_age
}
## Calculate time-specific weighted-mean for the IFR, weighted by current age-profile of infections
total_inf_by_age <- apply(inf_by_age, c("iteration", "date"), sum)
ifr.all <- array(apply(inf_by_age * ifr.tab, "age", `/`, total_inf_by_age), dim = dim(ifr.tab)[c("iteration", "date", "age")])
dimnames(ifr.all) <- list(iteration = dimnames(ifr.tab)$iteration, date = dimnames(ifr.tab)$date, age = dimnames(ifr.tab)$age)
ifr.eng <- apply(ifr.all, c("iteration", "date"), sum)
rm(inf_by_age, total_inf_by_age, ifr.all)

## Join England to the rest of the IFR array.
require(abind)
name.order <- c("iteration", "date", "age")
ifr.tab <- aperm(ifr.tab, name.order)
ifr.tab <- abind(ifr.eng, ifr.tab, along = 3); rm(ifr.eng)
names(dimnames(ifr.tab)) <- name.order
dimnames(ifr.tab)$age[1] <- "All ages"

require(cubelyr)
ifr.tab <- ifr.tab |>
    as.tbl_cube(met_name = "value") |>
    as_tibble() |>
    mutate(date = lubridate::as_date(date)) |>
    filter(date <= lubridate::ymd(date.data)) |>
    group_by(date, age) |>
    summarise(value = quantile(value, probs = c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) |>
    pivot_wider(id_cols = c(date, age), names_from = q, values_from = value)

write_csv(ifr.tab, paste0("IHR_", date.data, ".csv"))