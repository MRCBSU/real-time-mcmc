require(readr)
require(lubridate)

#######################################################################
## INPUT SETTINGS
#######################################################################

if(gp.flag){
	stop('No GP data for Italy')
} else {
    start.gp <- 1
    end.gp <- 1
}

## The 'hosp' stream in the code is linked to death data
if(!exists("hosp.flag")) hosp.flag <- 1	# 0 = off, 1 = on
if(hosp.flag){
    start.hosp <- 1
    ## Total days of data, or NULL to infer from length of file
    end.hosp <- NULL
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
first.matrix <- ymd(20200203)
day.num.first.matrix <- first.matrix - start.date + 1
cm.breaks <- seq(from = day.num.first.matrix, by = 7, length = 15)				# Day numbers where breaks happen
mat.dates <- start.date + cm.breaks - 1
lst <- readRDS(file.path(matrix.dir, "base_matrices.rds"))
lst$Lombardy$all$m <- lst$Lombardy$all$m * 1e7
cm.files <- "lombardy_8ag_contact.txt"
for(i in 1:length(cm.breaks))
    cm.files <- c(cm.files, paste0("lombardy_8ag_contact_ldwk", i, "_", google.data.date, ".txt"))
cm.bases <- file.path(proj.dir, "contact_mats", cm.files) ## Base matrices
cm.lockdown.fl <- paste0("Lombardy", mat.dates, "all.csv")
cm.lockdown <- file.path(matrix.dir, cm.lockdown.fl)
idx <- 1
if(!all(file.exists(cm.bases))){
    adf <- as.data.frame(lst$Lombardy$all$m)
    write_tsv(adf, cm.bases[idx], col_names = FALSE)
    for(fl in cm.lockdown){
        idx <- idx + 1
        mat <- read_csv(fl) * adf 
        write_tsv(mat, cm.bases[idx], col_names = FALSE)
    }
}
## Modifiers (which element of contact_parameters to use)
cm.mults <- file.path(proj.dir, "contact_mats", 
                      paste0("ag", nA, "_mult", c(0:1), ".txt"))
mult.order <- c(0)
for (d in mat.dates) {
	num.breakpoints.passed <- length(which(d >= m.param.breakpoints))
	mult.order <- c(mult.order, num.breakpoints.passed)
}
stopifnot(length(mult.order) == length(mat.dates) + 1)
if(!all(file.exists(cm.mults))){
    mult.mat <- lapply(unique(mult.order), function(x) matrix(x, nA, nA))
    for(i in 1:length(mult.mat)) write_tsv(as.data.frame(mult.mat[[i]]),
                                           cm.mults[i],
                                           col_names = FALSE)
    }
cm.mults <- cm.mults[mult.order+1]

## MCMC settings
num.iterations <- 200000
burnin <- 20000
adaptive.phase <- burnin / 2
thin.outputs <- 40 	# After how many iterations to output each set of NNI, deaths etc.
thin.params <- 20  # After how many iterations to output each set of parameters
stopifnot(thin.outputs %% thin.params == 0) # Need parameters on iterations we have outputs



############ NOTHING BELOW THIS LINE SHOULD NEED AMENDING WITH ANY REGULARITY ############
dir.data <- file.path(proj.dir, "data")
if(sys.nframe() <= 4){ ## Check if below source files might have already been loaded in
    source(file.path(proj.dir, "R/data/utils.R"))
    source(file.path(proj.dir, "config.R"))
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
hosp.data <- "NULL"
if (hosp.flag == 1) {
    hosp.data <- data.files
    if(!all(file.exists(hosp.data))) {
		print(hosp.data[which(!file.exists(hosp.data))])
        stop("Above hospitalisation data files does not exist")
	}
    if(is.null(end.hosp)) end.hosp <- set.end.date(end.hosp, hosp.data)
}
