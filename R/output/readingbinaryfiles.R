if(!exists("proj.dir")){
  proj.dir <- "/project/pandemic_flu/"
## proj.dir <- "/Volumes/Pandemic_flu/"
## proj.dir <- "~/bsu_pandemic/"
}

if(!exists("target.dir"))
  target.dir <- "./"

## output type
SMC.output <- FALSE

## run details - some default values
if(!exists("a")) a <- 7 ## number of age classes
if(!exists("d")) d <- 245 ## length in days of the run
if(!exists("i.saved")) i.saved <- 10000 ## number of iterations saved in coda files
if(!exists("i.summary")) i.summary <- 100 ## number of iterations of summary statistics stored on file
## regions <- c("London", "WestMidlands", "North", "South")
if(!exists("regions")) regions <- "London"
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
  if(SMC.output) states <- state.files <- vector("list", r)
  
  ## NNI files
  for(intr in 1:r)
    {
      NNI.files[[intr]] <- file(file.path(target.dir, paste0("NNI_", regions[intr])), "rb")
      if(SMC.output)
        ## state files
        state.files[[intr]] <- file(file.path(target.dir, paste0("state_", regions[intr])), "rb")
    }
  names(NNI.files) <- regions
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
      NNI[[intr]] <- readBin(NNI.files[[intr]], double(), n = i.summary * d * a)
      NNI[[intr]] <- array(NNI[[intr]], dim = c(a, d, i.summary))
    }
  names(NNI) <- regions
  
  lfx <- readBin(lfx.files, double(), n = i.saved)
  
  if(SMC.output)
    {
      for(intr in 1:r)
        states[[intr]] <- readBin(state.files[[intr]], double(), n = a * i.summary * 6) ## one each for S, E_1. E_2, I_1, I_2, p_lambda
      names(states) <- regions
    }
  
  for(inti in 1:npar){
	seek(coda.files[[inti]], -8, origin="end")
	num.params <- readBin(coda.files[[inti]], "integer")
	num.samples <- readBin(coda.files[[inti]], "integer")
	if (num.params != parameter.dims[inti]) {
		print(paste("WARNING: Expected", parameter.dims[inti], "chain(s) for", parameter.names[inti],
					"but", num.params, "found. Discarding the extra ones."))
	}
	expected.num.samples <- parameter.dims[inti] * i.saved
	if (num.samples != expected.num.samples) {
		print(paste("WARNING: Expected", expected.num.samples, "samples for", parameter.names[inti],
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
  for(intr in 1:r)
    close(NNI.files[[intr]])
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
