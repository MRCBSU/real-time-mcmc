require(readr)
require(lubridate)

#######################################################################
## INPUT SETTINGS
#######################################################################
if(gp.flag) {
    start.gp <- ll.start.date - start.date + 1			## What day to start running the likelihood on
    end.gp <- lubridate::as_date(date.data) - ll.reporting.delay - start.date + 1 ## Total days of data, or NULL to infer from length of file
} else {
    start.gp <- 1
    end.gp <- 1
}

## The 'hosp' stream in the code is linked to death data
if(!exists("hosp.flag") & !exists("adm.flag")) hosp.flag <- deaths.flag <- adm.flag <- 1 ## 0 = use hospital admissions, 1 = deaths
if(deaths.flag){
    start.hosp <- ifelse(data.desc == "reports", 35, 1) ## 35 # Day number on which to start likelihood calculation
    ## Total days of data, or NULL to infer from length of file
    end.hosp <- ifelse(use_deaths_up_to_now_flag, lubridate::as_date(date.data) - reporting.delay - start.date + 1, custom_deaths_end_date - start.date + 1)
    print(end.hosp)
    if(grepl("adjusted", data.desc)) end.hosp <- end.hosp + date.adj.data - lubridate::as_date(date.data)
} else if(adm.flag) {
    start.hosp <- 1
    end.hosp <- date.adm.str - start.date + 1
    if(cutoff_hosps_early & !deaths.flag & !hosp.flag) {
        end.hosp<- date_early_cutoff_hosps - start.date + 1
    }
} else {
    start.hosp <- 1
    end.hosp <- 1
}

## The 'sero' stream in the code
if(!exists("sero.flag")) sero.flag <- 1
## if(sero.flag){ ## Need to remove dependency  on rtm.plot as it may not necessarily be defined.
## 	if(exists("rtm.plot")) {
## 		start.sero <- min(rtm.plot$date) - start.date + 1
## 		end.sero <- max(rtm.plot$date) - start.date + 1
## 	} else {
## 		warning('Running sero likelihood between 14 and 85 days')
## 		start.sero <- 14
## 		end.sero <- 85 
## 	}
## } else {
## 	start.sero <- end.sero <- 1
## }

## The 'viro' stream in the code
viro.data <- NULL
viro.denom <- NULL

## The 'prev' stream in the code
if(!exists("prev.flag")) prev.flag <- 0
if(prev.flag){
    start.prev <- min(prev.lik.days)
    end.prev <- max(prev.lik.days)
} else {
    start.prev <- end.prev <- 1
}
# Vector of age-group descriptions
if(!exists("age.labs"))
    age.labs <- "All"


## CONTACT MATRICES SETTINGS
## Load Edwin's base matrices from contactsr
google.data.date_and_suff.str <- paste0(google.data.date, matrix.suffix) ## This line is confusing why would the date variable contain the suffix? suggest rename
matrix.dir <- file.path(
	proj.dir, "contact_mats",
	paste0("google_mobility_relative_matrices_", google.data.date_and_suff.str)
)
if (nA == 1) {
	cm.files <- rep("single_age.txt", length(cm.breaks) + 1)
	cm.bases <- file.path(proj.dir, "contact_mats", cm.files) ## Base matrices
} else {
  mat.dates <- start.date + cm.breaks - 1
  lst <- readRDS(file.path(matrix.dir, "base_matrices.rds"))
  lst$England$all$m <- lst$England$all$m * 1e7
  cm.files <- paste0("england_8ag_contact", google.data.date_and_suff.str, ".txt")
  for(i in 1:length(cm.breaks))
      cm.files <- c(cm.files, paste0("england_8ag_contact_ldwk", i, "_", google.data.date_and_suff.str, ifelse(flag.earlier_cm, "_earliercm_", ""), ".txt"))
  cm.bases <- file.path(proj.dir, "contact_mats", cm.files) ## Base matrices
  cm.lockdown.fl <- paste0("England", mat.dates, "all.csv")
  cm.lockdown <- file.path(matrix.dir, cm.lockdown.fl)
  idx <- 1
  if(!all(file.exists(cm.bases))){
      adf <- as.data.frame(lst$England$all$m)
      write_tsv(adf, cm.bases[idx], col_names = FALSE)
      for(fl in cm.lockdown){
          idx <- idx + 1
          mat <- read_csv(fl) * adf 
          write_tsv(mat, cm.bases[idx], col_names = FALSE)
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
} else if(contact.model == 4){
    cm.mults <- file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_mod4levels", 0:9, ".txt"))
    mult.order <- c(0, rep(1, length(cm.breaks)))
    mult.mat <- lapply(unique(mult.order), function(x){
        y <- (3*x)+(0:2)
        matrix(c(rep(y[2], 3 * nA), ## kids
                 rep(y[1], (nA - 4) * nA), ## adults except the very elderly
                 rep(y[3], nA)), nA, nA, byrow = TRUE)
        })
} else if(contact.model == 5){
    cm.mults <- file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_mod5levels", 0:9, ".txt"))
    mult.order <- c(0, rep(1, length(cm.breaks)))
    mult.mat <- lapply(unique(mult.order), function(x){
        y <- (4*x)+(0:3)
        matrix(c(rep(y[2], 3 * nA), ## kids,
                 rep(y[3], nA), ## university-aged
                 rep(y[1], (nA - 5) * nA), ## adults except the very elderly
                 rep(y[4], nA)), nA, nA, byrow = TRUE)
        })
} else if(contact.model == 6){ ## Each age group has a unique susceptibility
    cm.mults <- file.path(proj.dir, "contact_mats", paste0("ag", nA, "_mult_modAllLevels", 0:9, ".txt"))
    # If using modified contact matrices prior to the lockdown then use modified list to define cm_mults
    if(flag.earlier_cm) cm.breaks.strat.after.ld <- ifelse(cm.breaks > nday_lockdown_if_earlier_cm, 1, 0)
    mult.order <- c(0, `if`(flag.earlier_cm, cm.breaks.strat.after.ld, rep(1, length(cm.breaks))))
    mult.mat <- lapply(unique(mult.order), function(x){
        y <- ((nA-2)*x) + (0:(nA - 3))
        matrix(c(rep(y[2], 3 * nA), ## kids
                 rep(y[3], nA), ## 15-24
                 rep(y[1], nA), ## 25-44 (base age-group)
                 rep(y[4], nA), ## 45-64
                 rep(y[5], nA), ## 65-74
                 rep(y[6], nA)), ## 75+
               nA, nA, byrow = TRUE)
        })
}

if(!all(file.exists(cm.mults)))
    for(i in 1:length(mult.mat)) write_tsv(as.data.frame(mult.mat[[i]]),
                                       cm.mults[i],
                                       col_names = FALSE)
cm.mults <- cm.mults[mult.order+1]

## MCMC settings
num.iterations <- 5.5e6L
adaptive.phase <- 2e6L
burnin <- 4e6L
thin.outputs <- 600L ## After how many iterations to output each set of NNI, deaths etc.
thin.params <- 300L ## After how many iterations to output each set of parameters
# num.iterations <- 1e6L
# burnin <- 5e5L
# adaptive.phase <- 5e5L
# thin.outputs <- 500L ## After how many iterations to output each set of NNI, deaths etc.
# thin.params <- 250L ## After how many iterations to output each set of parameters
stopifnot(thin.outputs %% thin.params == 0) # Need parameters on iterations we have outputs
stored.covar <- 0
global.per.iter <- 1L



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
if (deaths.flag == 1) {
    hosp.data <- data.files
}
if (adm.flag) {
    hosp.data <- admsam.files
}
if(deaths.flag | adm.flag){
    if(!all(file.exists(hosp.data))) {
		print(hosp.data[which(!file.exists(hosp.data))])
        stop("Above hospitalisation data files does not exist")
	}
    if(is.null(end.hosp)) end.hosp <- set.end.date(end.hosp, hosp.data)
}

sero.data <- list(sample = "NULL",
                  positive = "NULL")

if (sero.flag == 1) {
    sero.data <- list(sample = serosam.files, positive = seropos.files)
    if(!all(sapply(sero.data, function(x) all(file.exists(x)))))
        stop("One of the specified serology data files does not exist")
    if(!format.inputs) { ## Need to define the limits for the serology likelihood
        ## Read in the samples file
        sero.lims <- do.call(bind_rows, lapply(serosam.files, read_tsv, col_names = FALSE)) %>%
            pivot_longer(cols = -1, names_to = "counts") %>%
            filter(value > 0)
        ## Find the date of the earliest and latest samples.
        start.sero <- min(sero.lims$X1) - start.date + 1
        end.sero <- ifelse(sero_cutoff_flag, sero.end.date - start.date + 1, max(sero.lims$X1) - start.date + 1)
    } else if(exists("rtm.plot")) {
        start.sero <- min(rtm.plot$date) - start.date + 1
        end.sero <- ifelse(sero_cutoff_flag, sero.end.date - start.date + 1, max(rtm.plot$date) - start.date + 1)
    } else if(exists("rtm.sam")) {
        start.sero <- min(rtm.sam$date) - start.date + 1
        end.sero <- ifelse(sero_cutoff_flag, sero.end.date - start.date + 1, max(rtm.sam$date) - start.date + 1)
    } else {
        warning('Running sero likelihood from day 1 to end\n')
        start.sero <- 1
        end.sero <- set.end.date(end.sero, sero.data)
    }
} else {
    start.sero <- end.sero <- 1
}

prev.data <- list(lmeans = "NULL",
                  lsds = "NULL")
if(prev.flag == 1){
    prev.data <- list(lmeans = prev.mean.files, lsds = prev.sd.files)
    if(!all(sapply(prev.data, function(x) all(file.exists(x)))))
        stop("One of the specified prevalence data files does not exist")
    if(is.null(end.prev)) end.prev <- set.end.date(end.prev, prev.data)
}
## Contact Model
if(!exists("cm.breaks")) {cm.breaks <- c(9, 16, 58, 72, 107, 114, 163, 212, 261, 268, 317)
    cm.bases <- file.path(proj.dir, "contact_mats", cm.bases)
    cm.mults <- file.path(proj.dir, "contact_mats", cm.mults)
}

num.threads <- nr * threads.per.regions

if (grepl("adjusted", data.desc)) {
	study_region_str <- "regions_hosp_aggregation = 5, 6, 7;"
} else {
	study_region_str <- ""
}

## DO WE WANT MCMC-STYLE CHAINS (0), OR SMC-STYLE PARTICLES (1)
mcmc.outs <- 0

# if(vacc.flag == 1) {
#     if(vac.n_doses == 2) {
#                     if(!all(file.exists(vac1.files)) || !all(file.exists(vacn.files))) {
#                 stop("One of the specified vaccination data files does not exist")
#             }
#     } else if(vac.n_doses == 3) {
#             if(!all(file.exists(vac1.files)) || !all(file.exists(vac2.files)) || !all(file.exists(vac3.files))) {
#                 stop("One of the specified vaccination data files does not exist")
#             }
#     }
# }
