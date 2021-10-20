prev.env <- new.env()
load.from <- file.path(previous.loc, "tmp.RData")
load(load.from, env = prev.env)
load.from <- file.path(previous.loc, "endstate.RData")
load(load.from, env = prev.env)
prev.params <- prev.env$params

value.eta.h <- prev.params$hosp_negbin_overdispersion[iteration.number.to.start.from,]
## value.dl <- latent period value is fixed
value.dI <- prev.params$infectious_period[iteration.number.to.start.from,]
if (prev.flag && prev.env$prev.flag) value.r1 <- prev.params$r1_period[iteration.number.to.start.from,]
## value.pgp GP stream not currently used
contact.reduction <- prev.params$contact_parameters[iteration.number.to.start.from,]
## if(contact.model == 6) { contact.reduction <- rep(contact.reduction, 2); contact.reduction[zero.contact.elements] <- 0 }
beta.rw.vals <- prev.params$log_beta_rw[iteration.number.to.start.from,]
beta.rw.vals <- add.extra.vals.per.region(beta.rw.vals, 0, nbetas)
if(nrow(beta.rw.vals) > nbetas)
    beta.rw.vals <- beta.rw.vals[c(1, 1 + sort(sample.int(nrow(beta.rw.vals)-1, nbetas-1))), ]
beta.rw.props <- add.extra.vals.per.region(prev.env$beta.rw.props, 0.02, nbetas)
if(nrow(beta.rw.props) > nbetas)
    beta.rw.props <- beta.rw.props[c(1, 1 + sort(sample.int(nrow(beta.rw.props)-1, nbetas-1))), ]
beta.rw.sd <- prev.params$log_beta_rw_sd[iteration.number.to.start.from, 2]
value.egr <- prev.params$exponential_growth_rate[iteration.number.to.start.from,]
value.nu <- prev.params$log_p_lambda_0[iteration.number.to.start.from,]
if(length(value.ifr) > ncol(prev.params$prop_case_to_hosp)){
    value.ifr <- c(prev.params$prop_case_to_hosp[iteration.number.to.start.from, ],
                   value.ifr[(ncol(prev.params$prop_case_to_hosp) + 1):length(value.ifr)])
} else value.ifr <- prev.params$prop_case_to_hosp[iteration.number.to.start.from, ]
## pars.dow Ignoring for now
if (!fix.sero.test.spec.sens && !prev.env$fix.sero.test.spec.sens) {
  sero.sens[1:ncol(prev.params$sero_test_sensitivity)] <- prev.params$sero_test_sensitivity[iteration.number.to.start.from,]
  sero.spec[1:ncol(prev.params$sero_test_sensitivity)] <- prev.params$sero_test_specificity[iteration.number.to.start.from,]
}
rm(prev.params)
rm(prev.env)
