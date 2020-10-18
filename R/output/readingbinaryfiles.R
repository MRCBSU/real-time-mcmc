if(!exists("proj.dir")){
  proj.dir <- "/project/pandemic_flu/"
## proj.dir <- "/Volumes/Pandemic_flu/"
## proj.dir <- "~/bsu_pandemic/"
}

if(!exists("target.dir"))
  target.dir <- out.dir

## output type
SMC.output <- FALSE

## run details - some default values
if(!exists("nA")) nA <- 1 ## number of age classes
if(!exists("ndays")) ndays <- 81 ## length in days of the run
if(!exists("i.saved")) i.saved <- 10000 ## number of iterations saved in coda files
if(!exists("i.summary")) i.summary <- 100 ## number of iterations of summary statistics stored on file
## regions <- c("London", "WestMidlands", "North", "South")
if(!exists("regions")) regions <- c("London", "Outside_London")
regions <- regions[!is.na(regions)]
r <- length(regions)
## ### ####### # #### ####### ######

## HAS THE PRIOR MODEL BEEN SPECIFIED? IF NO, NOTHING TO DO
if(exists("var.priors")){
  names(var.priors$distribution) <- names(var.priors$parameters) <- var.names
  
  
  stochastic.flags <- file.exists(file.path(target.dir, paste0("coda_", var.names)))
  
  parameter.names <- var.names[stochastic.flags]
  parameter.files <- file.path(target.dir, paste0("coda_", parameter.names))
  parameter.dims <- sapply(var.priors$distribution[stochastic.flags], length)
  
  npar <- length(parameter.files)
  
  ## ## OPEN FILE CONNECTIONS ## ##
  NNI <- NNI.files <- vector("list", r)
  dths.flag <- FALSE
  if(hosp.flag & !SMC.output) {
      Deaths <- Deaths.files <- vector("list", r)
      dths.flag <- !dths.flag
  }
  cases.flag <- FALSE
  if(gp.flag & !SMC.output) {
      Cases <- Cases.files <- vector("list", r)
      cases.flag <- !cases.flag
  }
  if(!exists("prev.flag")) prev.flag <- FALSE
  if(prev.flag & !SMC.output)
      Prev <- Prev.files <- vector("list", r)
  if(SMC.output) states <- state.files <- vector("list", r)
  
  ## NNI files
  for(intr in 1:r)
    {
        NNI.files[[intr]] <- file(file.path(target.dir, paste0("NNI_", regions[intr])), "rb")
        if(dths.flag)
            Deaths.files[[intr]] <- file(file.path(target.dir, paste0("Hosp_", regions[intr])), "rb")
        if(cases.flag)
            Cases.files[[intr]] <- file(file.path(target.dir, paste0("GP_", regions[intr])), "rb")
        if(prev.flag)
            Prev.files[[intr]] <- file(file.path(target.dir, paste0("Prev_", regions[intr])), "rb")
        if(SMC.output)
            ## state files
            state.files[[intr]] <- file(file.path(target.dir, paste0("state_", regions[intr])), "rb")
    }
  names(NNI.files) <- regions
  if(dths.flag) names(Deaths.files) <- regions
  if(cases.flag) names(Cases.files) <- regions
  if(prev.flag) names(Prev.files) <- regions
  if(SMC.output) names(state.files) <- regions
  
  ## lfx files
  lfx.files <- file(file.path(target.dir, "coda_lfx"), "rb")
  ## coda files
  params <- coda.files <- list()
  for(inti in 1:npar)
    coda.files[[inti]] <- file(parameter.files[inti], "rb")
  ## ## ## ##
  
  ## ## READ IN FROM FILE CONNECTIONS ## ##
  for(intr in 1:r)
    {
      NNI[[intr]] <- readBin(NNI.files[[intr]], double(), n = i.summary * ndays * nA)
      NNI[[intr]] <- array(NNI[[intr]], dim = c(nA, ndays, i.summary))
      if(dths.flag){
          Deaths[[intr]] <- readBin(Deaths.files[[intr]], double(), n = i.summary * ndays * nA) %>%
              array(dim = c(nA, ndays, i.summary))
      }
      if(cases.flag){
          Cases[[intr]] <- readBin(Cases.files[[intr]], double(), n = i.summary * ndays * nA) %>%
              array(dim = c(nA, ndays, i.summary))
      }
      if(prev.flag){
          Prev[[intr]] <- readBin(Prev.files[[intr]], double(), n = i.summary * ndays * nA) %>%
              array(dim = c(nA, ndays, i.summary))
      }
    }
  names(NNI) <- regions
  if(dths.flag) names(Deaths) <- regions
  if(cases.flag) names(Cases) <- regions
  if(prev.flag) names(Prev) <- regions
  
  lfx <- readBin(lfx.files, double(), n = i.saved)
  
  if(SMC.output)
  {
      for(intr in 1:r)
          states[[intr]] <- readBin(state.files[[intr]], double(), n = nA * i.summary * 6) ## one each for S, E_1. E_2, I_1, I_2, p_lambda
      names(states) <- regions
    }
  
  for(inti in 1:npar){
	seek(coda.files[[inti]], -8, origin="end")
	num.params <- readBin(coda.files[[inti]], "integer")
	num.samples <- readBin(coda.files[[inti]], "integer")
	if (num.params != parameter.dims[inti]) {
		print(paste("WARNING: Expected", parameter.dims[inti], "chain(s) for", parameter.names[inti],
                            "but", num.params, "found. Discarding extra ones and not plotting missing ones."))
	}
	## expected.num.samples <- parameter.dims[inti] * i.saved
	if (num.samples != i.saved) {
		print(paste("WARNING: Expected", i.saved, "samples for", parameter.names[inti],
					"but", num.samples, "found. Discarding the extra ones."))
	}
	seek(coda.files[[inti]], 0)
	params[[inti]] <- readBin(coda.files[[inti]], "double", num.samples * num.params)
 	params[[inti]] <- array(params[[inti]], dim = c(num.params, num.samples))
	params[[inti]] <- params[[inti]][1:parameter.dims[inti], 1:i.saved, drop=FALSE]
  }
  names(params) <- parameter.names
  ## ## ## ##
  
  ## ## CLOSE FILE CONNECTIONS ## ##
  for(intr in 1:r){
      close(NNI.files[[intr]])
      if(dths.flag) close(Deaths.files[[intr]])
      if(cases.flag) close(Cases.files[[intr]])
      if(prev.flag) close(Prev.files[[intr]])
  }
  close(lfx.files)
  if(SMC.output)
    for(intr in 1:r)
      close(state.files[[intr]])
  for(inti in 1:npar)
    close(coda.files[[inti]])
  ## ## ## ## ## ## ## ## ## ## ## #
  
} else {
  
  cat("Nothing done: no prior model specified\n")
  
}
