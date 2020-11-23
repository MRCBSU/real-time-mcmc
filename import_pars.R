add.extra.vals.per.region <- function(vec, val, num) {
  mat <- matrix(vec, ncol = nr)
  rows.to.add <- nbetas - length(vec) / nr 
  mat.new <- matrix(val, nrow = rows.to.add, ncol = nr)
  return(rbind(mat, mat.new))
}
  

prev.env <- new.env()
load.from <- file.path(previous.run.to.use, "tmp.RData")
load(load.from, env = prev.env)
load.from <- file.path(previous.run.to.use, "mcmc.RData")
load(load.from, env = prev.env)
prev.params <- prev.env$params

value.eta.h <- prev.params$hosp_negbin_overdispersion[iteration.number.to.start.from,]
#value.dl <- latent period value is fixed
value.dI <- prev.params$infectious_period[iteration.number.to.start.from,]
if (prev.flag && prev.env$prev.flag) value.dR <- prev.params$r1_period[iteration.number.to.start.from,]
#value.pgp GP stream not currently used
contact.reduction <- prev.params$contact_parameters[iteration.number.to.start.from,]
beta.rw.vals <- prev.params$log_beta_rw[iteration.number.to.start.from,]
beta.rw.vals <- add.extra.vals.per.region(beta.rw.vals, 0, nbetas * nr)
beta.rw.sd <- prev.params$log_beta_rw_sd[iteration.number.to.start.from, 2]
value.egr <- prev.params$exponential_growth_rate[iteration.number.to.start.from,]
value.nu <- prev.params$log_p_lambda_0[iteration.number.to.start.from,]
value.ifr <- prev.params$prop_case_to_hosp[iteration.number.to.start.from,]
#pars.dow Ignoring for now
if (!fix.sero.test.spec.sens && !prev.env$fix.sero.test.spec.sens) {
  sero.sens <- prev.params$sero_test_sensitivity[iteration.number.to.start.from,]
  sero.spec <- prev.params$sero_test_specificity[iteration.number.to.start.from,]
}

