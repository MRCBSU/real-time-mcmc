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
r <- length(regions)
## ### ####### # #### ####### ######

## HAS THE PRIOR MODEL BEEN SPECIFIED? IF NO, NOTHING TO DO
if(exists("var.priors")){
  names(var.priors$distribution) <- names(var.priors$parameters) <- var.names
  
  
  stochastic.flags <- file.exists(paste(target.dir, "coda_", var.names, sep = ""))
  
  parameter.names <- var.names[stochastic.flags]
  parameter.files <- paste(target.dir, "coda_", parameter.names, sep = "")
  parameter.dims <- sapply(var.priors$distribution[stochastic.flags], length)
  
  npar <- length(parameter.files)
  
  ## ## OPEN FILE CONNECTIONS ## ##
  NNI <- NNI.files <- vector("list", r)
  if(SMC.output) states <- state.files <- vector("list", r)
  
  ## NNI files
  for(intr in 1:r)
    {
      NNI.files[[intr]] <- file(paste(target.dir, "NNI_", regions[intr], sep = ""), "rb")
      if(SMC.output)
        ## state files
        state.files[[intr]] <- file(paste(target.dir, "state_", regions[intr], sep = ""), "rb")
    }
  names(NNI.files) <- regions
  if(SMC.output) names(state.files) <- regions
  
  ## lfx files
  lfx.files <- file(paste(target.dir, "coda_lfx", sep = ""), "rb")
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
    params[[inti]] <- readBin(coda.files[[inti]], double(), n = i.saved * parameter.dims[inti])
    params[[inti]] <- array(params[[inti]], dim = c(parameter.dims[inti], i.saved))
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
