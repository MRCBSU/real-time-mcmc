drw <- function(iter, pars, nc = 1){

    ## Extract the iteration
    parsx <- lapply(pars, function(x) x[iter, ])

    ## Extract the parameters required to evaluate the random-walk density
    beta <- parsx$log_beta_rw
    sigma <- parsx$log_beta_rw_sd[nc]
    beta <- beta[beta != 0]

    sum(dnorm(beta, sd = sigma, log = TRUE))

}
