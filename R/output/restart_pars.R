require(rmarkdown)
require(knitr)

load("tmp.RData")
load("mcmc.RData")

if("hosp_negbin_overdispersion" %in% names(params))
    value.eta.h <- params$hosp_negbin_overdispersion[i.saved,]
## value.dl <- latent period value is fixed
if("infectious_period" %in% names(params))
    value.dI <- params$infectious_period[i.saved,]
if (prev.flag)
    value.r1 <- params$r1_period[i.saved,]
## value.pgp GP stream not currently used
if("contact_parameters" %in% names(params))
    contact.reduction <- params$contact_parameters[i.saved,]
## if(contact.model == 6) { contact.reduction <- rep(contact.reduction, 2); contact.reduction[zero.contact.elements] <- 0 }
if("log_beta_rw" %in% names(params)){
    beta.rw.vals <- params$log_beta_rw[i.saved,]
    beta.rw.vals <- add.extra.vals.per.region(beta.rw.vals, 0, nbetas)
    if(nrow(beta.rw.vals) > nbetas)
        beta.rw.vals <- beta.rw.vals[c(1, 1 + sort(sample.int(nrow(beta.rw.vals)-1, nbetas-1))), ]
    beta.rw.props <- add.extra.vals.per.region(beta.rw.props, 0.02, nbetas)
    if(nrow(beta.rw.props) > nbetas)
        beta.rw.props <- beta.rw.props[c(1, 1 + sort(sample.int(nrow(beta.rw.props)-1, nbetas-1))), ]
}
if("log_beta_rw_sd" %in% names(params))
    beta.rw.sd <- params$log_beta_rw_sd[i.saved, 2]
if("exponential_growth_rate" %in% names(params))
    value.egr <- params$exponential_growth_rate[i.saved,]
if("log_p_lambda_0" %in% names(params))
    value.nu <- params$log_p_lambda_0[i.saved,]
if(length(value.ifr) > ncol(params$prop_case_to_hosp)){
    value.ifr <- c(params$prop_case_to_hosp[i.saved, ],
                   value.ifr[(ncol(params$prop_case_to_hosp) + 1):length(value.ifr)])
} else value.ifr <- params$prop_case_to_hosp[i.saved, ]
## pars.dow Ignoring for now
if ("sero_test_sensitivity" %in% names(params)){
  sero.sens[1:ncol(params$sero_test_sensitivity)] <- params$sero_test_sensitivity[i.saved,]
  sero.spec[1:ncol(params$sero_test_sensitivity)] <- params$sero_test_specificity[i.saved,]
}

pars.template.loc <- file.path(proj.dir, "inputs", "mod_pars.Rmd")
knit(input = pars.template.loc, output = file.path(out.dir, "mod_pars_new.txt"))
