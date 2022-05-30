## Written for inclusion in report-updated.Rmd

## Presumes tmp.RData mcmc.RData and output_matrices.RData are all loaded in to the workspace.
if(!exists("out.dir"))
    load("tmp.RData") ## There is a function somewhere in the repo that will guess sensible (relative) file paths. Update to use this.
if(!exists("ifr"))
    load(file.path(out.dir, "output_matrices.RData"))

## ## Base R proved quicker than tidyr equivalents, to...
## Aggregate infections over region
## Aggegrate the two youngest age groups
sum.up.and.merge <- function(inf){
		 inf_by_age <- apply(inf, c("age", "date", "iteration"), sum)
		 inf_by_age <- merge.youngest.age.groups(inf_by_age, idx = c("iteration", "date"))
		names(dimnames(inf_by_age))[1] <- "age"
		dimnames(inf_by_age)$age <- dimnames(ifr)$age
		inf_by_age
}


## Only look at IFR iterations for which the infections array has been evaluated
ifr.tab <- extract(ifr, indices=list(parameter.to.outputs), dims="iteration")

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
