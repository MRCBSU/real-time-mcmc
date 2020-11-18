prev.env <- new.env()
load.from <- file.path(previous.run.to.use, "tmp.RData")
load(load.from, env = prev.env)
load.from <- file.path(previous.run.to.use, "mcmc.RData")
load(load.from, env = prev.env)
prev.params <- prev.env$params

value.eta.h <- prev.params$hosp_negbin_overdispersion[iteration.number.to.start.from,]
#value.dl <- latent period value is fixed
value.dI <- prev.params$infectious_period[iteration.number.to.start.from,]
value.dR <- prev.params$r1_period[iteration.number.to.start.from,]
#value.pgp GP stream not currently used
contact.reduction <- prev.params$contact_parameters[iteration.number.to.start.from,]
beta.rw.vals <- prev.params$log_beta_rw[iteration.number.to.start.from,]
beta.rw.sd <- prev.params$log_beta_sd[iteration.number.to.start.from,]
value.egr <- prev.params$exponential_growth_rate[iteration.number.to.start.from,]
value.nu <- prev.params$log_p_lambda_0[iteration.number.to.start.from,]
value.ifr <- prev.params$prop_case_to_hosp[iteration.number.to.start.from,]
#pars.dow Ignoring for now
if (!fix.sero.test.spec.sens && !prev.env$fix.sero.test.spec.sens) {
  # TODO: fix this
  stop("Not implemented importing serology fields yet")
}
stopifnot(all(!is.na(beta.rw.vals)))
stopifnot(length(beta.rw.vals) == nbetas*nr)
stopifnot(all(beta.rw.vals[static.zero.beta.locs] == 0))
stopifnot(all(!is.na(beta.rw.props)))
stopifnot(length(beta.rw.props) == nbetas*nr)
stopifnot(all(beta.rw.props[static.zero.beta.locs] == 0))
stopifnot(all(beta.rw.props[-static.zero.beta.locs] > 0))

