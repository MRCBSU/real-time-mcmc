library(lubridate)
library(tidyverse)

#######################################################################
## INPUT SETTINGS
#######################################################################

if(gp.flag){
    start.gp <- 15			# What day to start running the likelihood on
    end.gp <- NULL			# Total days of data, or NULL to infer from length of file
} else {
    start.gp <- 1
    end.gp <- 1
}

## The 'hosp' stream in the code is linked to death data
if(!exists("hosp.flag")) hosp.flag <- 1	# 0 = off, 1 = on
if(hosp.flag){
    start.hosp <- ifelse(data.desc == "reports", 35, 1) ## 35 # Day number on which to start likelihood calculation
    ## Total days of data, or NULL to infer from length of file
    end.hosp <- lubridate::as_date(date.data) - reporting.delay - start.date + 1
} else {
    start.hosp <- 1
    end.hosp <- 1
}
## The 'sero' stream in the code
if(!exists("sero.flag")) sero.flag <- 1
if(sero.flag){ ## Need to remove dependency  on rtm.plot as it may not necessarily be defined.
	if(exists("rtm.plot")) {
		start.sero <- min(rtm.plot$date) - start.date + 1
		end.sero <- max(rtm.plot$date) - start.date + 1
	} else {
		warning('Running sero likelihood for whole period')
		start.sero <- 1
		end.sero <- ndays 
	}
} else {
	start.sero <- end.sero <- 1
}
## The 'viro' stream in the code
viro.data <- NULL
viro.denom <- NULL

# Vector of age-group descriptions
if(!exists("age.labs"))
    age.labs <- "All"


## CONTACT MATRICES SETTINGS
## Load Edwin's base matrices from contactsr
matrix.dir <- file.path(
	proj.dir, "contact_mats",
	paste0("google_mobility_relative_matrices_", google.data.date)
)
last.break <- as.integer(ymd(google.data.date) - ymd(start.date), unit = "days") - 4
cm.breaks <- seq(from = 36, to = last.break, by = 7)
lst <- readRDS(file.path(matrix.dir, "base_matrices.rds"))
if (running.England) {
  cm.region.name <- "England"
} else if (nr > 1) {
  stop("Only support multiple regions for England")
} else if (regions == "Wales") {
  cm.region.name <- "Wales"
} else if (regions == "Northern_Ireland") {
  cm.region.name <- "Northern Ireland"
} else if (regions == "Scotland") {
  cm.region.name <- "Scotland"
} else {
  stop("Unknown region")
}
adf <- as.data.frame(lst[[cm.region.name]]$all$m * 1e7)
if (nA == 1 && !include.google) {
	cm.files <- rep("single_age.txt", length(cm.breaks) + 1)
	cm.bases <- file.path(proj.dir, "contact_mats", cm.files) ## Base matrices
} else {
  mat.dates <- start.date + cm.breaks - 1
  cm.out.region <- str_replace_all(cm.region.name, " ", "_")
  cm.file.base <- paste0(cm.out.region, "_", nA, "ag_contact")
  cm.files <- paste0(cm.file.base, ".txt")
  cm.run.stamp <- google.data.date
  for(i in 1:length(cm.breaks)) {
      cm.files <- c(cm.files, paste0(cm.file.base, "_ldwk", i, "_", cm.run.stamp, ".txt"))
  }
  cm.bases <- file.path(proj.dir, "contact_mats", cm.files) ## Base matrices
  cm.lockdown.fl <- paste0(cm.region.name, mat.dates, "all.csv")
  cm.lockdown <- file.path(matrix.dir, cm.lockdown.fl)
  stopifnot(length(cm.bases) == length(cm.lockdown) + 1)
  if (create.counterfactual) {
    max.matrix.date <- intervention.date - 7
    mats.to.use <- mat.dates <= max.matrix.date
    mats.to.use <- c(TRUE, mats.to.use) # Always use first pre-lockdown matrix
    last.mat.to.use <- max(which(mats.to.use))
    cm.bases[!mats.to.use] <- cm.bases[last.mat.to.use]
  }
  idx <- 1
  if(!all(file.exists(cm.bases))){
	 # Pre-lockdown matrix (cm.bases[1])
     if (nA == 1) {
		R <- max(abs(eigen(adf, only.values = TRUE)$value))
		cat(R, file = cm.bases[idx])
	  } else {
		write_tsv(adf, cm.bases[idx], col_names = FALSE)
	  }
      # Post-lockdown matrices
      for(fl in cm.lockdown){
          idx <- idx + 1
          mat <- read_csv(fl) * adf 
		  if (nA == 1) {
	        R <- max(abs(eigen(mat, only.values = TRUE)$value))
		    cat(R, file = cm.bases[idx])
		  } else {
            write_tsv(mat, cm.bases[idx], col_names = FALSE)
		  }
      }
   }
}
## Modifiers (which element of contact_parameters to use)
if(contact.model == 1){
    cm.mults <- file.path(proj.dir, "contact_mats", 
                          paste0("ag", nA, "_mult", 0:9, ".txt"))
    mult.order <- c(0, rep(1, length(cm.breaks)))
    ## mult.order <- 0:length(cm.breaks)
    mult.mat <- lapply(unique(mult.order), function(x) matrix(x, nA, nA))
} else if(contact.model == 2){
    cm.mults <- file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_3levels", 0:9, ".txt"))
    mult.order <- c(0, rep(1, length(cm.breaks)))
    mult.mat <- lapply(unique(mult.order), function(x){
        y <- (2*x)-(1:0)
        if(x==0) y <- rep(0, 2)
        matrix(c(rep(y[1], nA * (nA - 1)),
               rep(y[2], nA)), nA, nA, byrow = TRUE)
    })
} else if(contact.model == 3){
    cm.mults <- file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_mod3levels", 0:9, ".txt"))
    mult.order <- c(0, rep(1, length(cm.breaks)))
    mult.mat <- lapply(unique(mult.order), function(x){
        y <- (2*x)+(0:1)
        matrix(c(rep(y[1], nA * (nA - 1)),
                 rep(y[2], nA)), nA, nA, byrow = TRUE)
    })
}
if(!all(file.exists(cm.mults)))
    for(i in 1:length(mult.mat)) write_tsv(as.data.frame(mult.mat[[i]]),
                                       cm.mults[i],
                                       col_names = FALSE)
cm.mults <- cm.mults[mult.order+1]

## MCMC settings
num.iterations <- 800e3
stopifnot(num.iterations < 1e6) # mod_inputs.txt format does not support integers >= one million
burnin <- 30e3
adaptive.phase <- burnin / 2
thin.outputs <- 30 	# After how many iterations to output each set of NNI, deaths etc.
thin.params <- 10  # After how many iterations to output each set of parameters
stopifnot(thin.outputs %% thin.params == 0) # Need parameters on iterations we have outputs



############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############
dir.data <- file.path(proj.dir, "data")
if(sys.nframe() <= 4){ ## Check if below source files might have already been loaded in
    source(file.path(proj.dir, "R/data/utils.R"))
    source(file.path(proj.dir, "config.R"))
}

## Map what we call regions (LHS) to the NHS region(s) they contain
## These no longer calculated `on the fly' and should be handled within the data/population folder.
## Use objects nhs.regions and pop
load(build.data.filepath("population", "pop_nhs.RData"))

get.nhs.region <- function(reg, rlist = nhs.regions){
    if(reg %in% names(nhs.regions)){
        return(reg)
    } else if(toupper(reg) %in% names(nhs.regions)) return(toupper(reg))
}
## Check that regions have population specified
for (region in regions) {
    if (!get.nhs.region(region) %in% names(nhs.regions)) {
        stop(paste(region, "is not specified in `nhs.regions`. Options are:",
                   paste0(names(nhs.regions), collapse=", ")))
    }
}

# If end.gp and/or end.hosp are none then read from data files
set.end.date <- function(user.value, data.file) {
	if (is.null(user.value)) {
		return(max(length(readLines(data.file[1]))))
	} else {
		return(user.value)
	}
}
# Where are the data files?
dir.data <- file.path(proj.dir, "data")
source(file.path(proj.dir, "R/data/utils.R"))

gp.data <- "NULL"
gp.denom <- "NULL"
if (gp.flag == 1) {
	gp.data <- build.data.filepath("RTM_format", "linelist", date.data, ".txt")
	gp.denom <- build.data.filepath("RTM_format", "ll_denom", date.data, ".txt")
	if(is.null(end.gp)) end.gp <- set.end.date(end.gp, gp.data)
}
hosp.data <- "NULL"
if (hosp.flag == 1) {
    hosp.data <- data.files
    if(!all(file.exists(hosp.data))) {
		print(hosp.data[which(!file.exists(hosp.data))])
        stop("Above hospitalisation data files does not exist")
	}
    if(is.null(end.hosp)) end.hosp <- set.end.date(end.hosp, hosp.data)
}
if (create.counterfactual) {
  max.hosp.date <- intervention.date + 7
  max.hosp.time <- max.hosp.date - start.date
  end.hosp <- min(end.hosp, max.hosp.time)
}
sero.data <- list(sample = "NULL",
                  positive = "NULL")
if (sero.flag == 1) {
    sero.data <- list(sample = serosam.files, positive = seropos.files)
    if(!all(sapply(sero.data, function(x) all(file.exists(x)))))
        stop("One of the specified serology data files does not exist")
    if(is.null(end.sero)) end.sero <- set.end.date(end.sero, sero.data)
}
## Contact Model
if(!exists("cm.breaks")) {cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
    cm.bases <- file.path(proj.dir, "contact_mats", cm.bases)
    cm.mults <- file.path(proj.dir, "contact_mats", cm.mults)
}

num.threads <- nr

if (data.desc == "adjusted") {
	study_region_str <- "regions_hosp_aggregation = 5, 6, 7;"
} else {
	study_region_str <- ""
}
