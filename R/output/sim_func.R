sim_rtm <- function(iter, rtm.exe = Sys.info()["nodename"]){

    ## Fix output location
    out.loc <- paste0(projections.basedir, "_", Sys.getpid())
    if(!file.exists(out.loc)){
        system(paste("cp -r",
                     projections.basedir,
                     out.loc))
    }
    fl.pars <- file.path(out.loc, "sim_mod_pars.Rmd")
    if(!file.exists(fl.pars))
        file.copy(gsub("mod_pars", "sim_mod_pars", pars.template.loc),
                  fl.pars)

    Deaths <- NNI <- Cases <- Prevs <- vector("list", nr)
    Deaths.files <- Prev.files <- NNI.files <- Cases.files <- vector("list", nr)
    
    ## Set-up output objects
    ## Get the parameters into the right variable names.
    if(gp.flag) value.eta <- params$gp_negbin_overdispersion[iter, , drop = FALSE]
    if(hosp.flag) value.eta.h <- params$hosp_negbin_overdispersion[iter, , drop = FALSE]
    value.dI <- params$infectious_period[iter, , drop = FALSE]
    if(gp.flag) value.pgp  <- params$prop_case_to_GP_consultation[iter, , drop = FALSE]
    contact.reduction <- params$contact_parameters[iter, , drop = FALSE]
    beta.rw.vals <- params$log_beta_rw[iter, , drop = FALSE]
    log_beta_rw_sd <- params$log_beta_rw_sd[iter, -1, drop = FALSE]
    value.egr <- params$exponential_growth_rate[iter, , drop = FALSE]
    value.nu <- params$log_p_lambda_0[iter, , drop = FALSE]
    if(hosp.flag) value.ifr <- params$prop_case_to_hosp[iter, , drop = FALSE]
    if(gp.flag) pars.dow <- params$day_of_week_effects[iter, , drop = FALSE]
    sero.sens <- params$sero_test_sensitivity[iter, , drop = FALSE]
    sero.spec <- params$sero_test_specificity[iter, , drop = FALSE]
    
    ## ## Compile the mod_pars.txt file
    ## render(fl.pars, output_file = "mod_pars.txt",
    ##        output_dir = out.loc, output_format = plain_document, quiet = TRUE)
    knit(input = fl.pars, output = file.path(out.loc, "mod_pars.txt"), quiet = TRUE)
    
    ## Run the code
    setwd(out.loc)
    exit_code <- system(file.path(proj.dir, paste0("rtm_", rtm.exe)), intern = FALSE)
    if(exit_code != 0) stop(paste("Error running in", out.loc, "on iteration", iter))
    
    ## Read the outputs in and append to output objects
    for(intr in 1:nr)
    {
        NNI.files[[intr]] <- file(paste0("NNI_", regions[intr]), "rb")
        if(dths.flag)
            Deaths.files[[intr]] <- file(paste0("Hosp_", regions[intr]), "rb")
        if(cases.flag)
            Cases.files[[intr]] <- file(paste0("GP_", regions[intr]), "rb")
        if(prev.flag)
            Prev.files[[intr]] <- file(paste0("Prev_", regions[intr]), "rb")
    }
    names(NNI.files) <- regions
    if(dths.flag) names(Deaths.files) <- regions
    if(cases.flag) names(Cases.files) <- regions
    if(prev.flag) names(Prev.files) <- regions
    
    for(intr in 1:r)
    {
        NNI[[intr]] <- readBin(NNI.files[[intr]], double(), n = num.iterations * ndays * nA) %>%
            array(dim = c(nA, ndays, num.iterations))
        close(NNI.files[[intr]])
        if(dths.flag){
            Deaths[[intr]] <- readBin(Deaths.files[[intr]], double(), n = num.iterations * ndays * nA) %>%
                      array(dim = c(nA, ndays, num.iterations))
            close(Deaths.files[[intr]])
        }
        if(cases.flag){
            Cases[[intr]] <- readBin(Cases.files[[intr]], double(), n = num.iterations * ndays * nA) %>%
                      array(dim = c(nA, ndays, num.iterations))
            close(Cases.files[[intr]])
        }
        if(prev.flag){
            Prevs[[intr]] <- readBin(Prev.files[[intr]], double(), n = num.iterations * ndays * nA) %>%
                array(dim = c(nA, ndays, num.iterations))
            close(Prev.files[[intr]])
        }
    }
    ## if((100 * iter / niter) > (pct + 1)){
    ##     pct <- pct + 1
    ##     cat(pct, "% Complete\n")
    ## }
    setwd("..")
    list(NNI = NNI, Deaths = Deaths, Cases = Cases, Prevs = Prevs)
}
