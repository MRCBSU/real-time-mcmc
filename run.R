library(rmarkdown)
library(tidyverse)

## Location of this script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

## Where are various directories?
file.loc <- dirname(thisFile())
#file.loc = "~/real-time-mcmc"
proj.dir <- file.loc
source(file.path(proj.dir, "config.R"))
source(file.path(proj.dir, "R/data/utils.R"))

system(paste("mkdir -p", out.dir))

## do we need to do formatting?
format.inputs <- TRUE

## Will code need to be recompiled?
compile.code <- FALSE

## Do we want to actually run the code?
run.code <- FALSE

## Do we want to automatically run post-processing R code?
run.outputs <- FALSE

## Which code is being considered
if(!exists("gp.flag")) gp.flag <- 1
if(!exists("hosp.flag")) hosp.flag <- deaths.flag <- 1
if(!exists("adm.flag")) adm.flag <- !hosp.flag ## Not typically set by config.R - check
if(!exists("sero.flag")) sero.flag <- 1
if(!exists("viro.flag")) viro.flag <- 0
if(!exists("prev.flag")) prev.flag <- 0
if(!exists("NHSBT.flag")) NHSBT.flag <- 1 # NHSBT == 1, RCGP == 0
if(!exists("RocheS.flag")) RocheS.flag <- 1 # RocheS == 1, RocheN == 0

if(!adm.flag) admsam.files <- NULL ## The output filenames for the data will be determined within format_hosp_admissions.R

if (region.type == "NHS") {
	source(file.path(proj.dir, "R/data/get_NHS_pop.R"))
} else if (region.type == "ONS") {
	source(file.path(proj.dir, "R/data/get_ONS_pop.R"))
} else stop("Unknown region type for population")


# Moved serology calls to the top of the run script
if(sero.flag){
    str.collect <- ifelse(NHSBT.flag, "NHSBT", "RCGP")
    serosam.files <- paste0(data.dirs["sero"], "/", sero.end.date, "_", regions, "_", nA, "ag_", str.collect, "samples", ifelse(!sero_cutoff_flag, "", paste0("_dropsero_", gsub("-", "",toString(sero.end.date)))), ".txt")
    seropos.files <- paste0(data.dirs["sero"], "/", sero.end.date, "_", regions, "_", nA, "ag_", str.collect, "positives", ifelse(!sero_cutoff_flag, "", paste0("_dropsero_", gsub("-", "",toString(sero.end.date)))), ".txt")
    if(!format.inputs) {
        # Copy files but take into account possible name differences due to custom cutoff
        file.copy(file.path(data.dirs["sero"], paste0("sero_samples_data", ifelse(!sero_cutoff_flag, "", paste0("_dropsero_", gsub("-", "",toString(sero.end.date)))), ".csv")), out.dir)
        file.rename(file.path(out.dir, paste0("sero_samples_data", ifelse(!sero_cutoff_flag, "", paste0("_dropsero_", gsub("-", "",toString(sero.end.date)))), ".csv")), file.path(out.dir, "sero_samples_data.csv"))
        file.copy(file.path(data.dirs["sero"], paste0("sero_positives_data", ifelse(!sero_cutoff_flag, "", paste0("_dropsero_", gsub("-", "",toString(sero.end.date)))), ".csv")), out.dir)
        file.rename(file.path(out.dir, paste0("sero_positives_data", ifelse(!sero_cutoff_flag, "", paste0("_dropsero_", gsub("-", "",toString(sero.end.date)))), ".csv")), file.path(out.dir, "sero_positives_data.csv"))
        names(serosam.files) <- names(seropos.files) <- regions
    }
} else {
  serosam.files <- seropos.files <- NULL
}

## If these files don't already exits, make them
if(deaths.flag){
    data.files <- paste0(data.dirs["deaths"], "/",
                         data.desc,
                         date.data, "_",
                         regions, "_",
                         nA, "ag",
                         ifelse(flg.confirmed, "CONF", ""),
                         reporting.delay, "delay",
                         # Take into account possible name differences due to custom cutoff
                         ifelse(use_deaths_up_to_now_flag, "", paste0("_", custom_deaths_end_date))
                         )
    if (exists("flg.cutoff")){
        if(flg.cutoff)
            data.files <- paste0(data.files, "cutoff", str.cutoff)
    }
    data.files <- paste0(data.files, ".txt")
    names(data.files) <- regions
}
if(gp.flag){
    cases.files <- paste0(data.dirs["cases"], "/", date.data, "_", regions, "_", nA, "_pillar_2_", ifelse(symptoms, "symptoms", "all"), ".txt")
    denoms.files <- paste0(data.dirs["cases"], "/", date.data, "_", regions, "_", nA, "_popdenom.txt")
} else {
    cases.files <- NULL
    denoms.files <- NULL
}
if(prev.flag){
#    prev.file.prefix <- paste0(data.dirs["prev"], "/date_prev", "_")
    prev.file.txt <- ifelse(all(diff(prev.lik.days) == 1),
                            #paste(min(prev.lik.days), "every_day", max(prev.lik.days)-300, sep = "_"),
                            paste(min(prev.lik.days), "every_day",456, sep = "_"),
                            paste0(prev.lik.days[c(1:which(prev.lik.days==456))], collapse = "_"))
 #                           paste0(prev.lik.days, collapse = "_"))
    if (exists("date.prev")) {
		prev.file.prefix <- paste0(data.dirs["prev"], "/", date.prev, "_", prev.file.txt, "_")
                          	}
        else {
		prev.file.prefix <- paste0(data.dirs["prev"], "/date_prev_", prev.file.txt, "_")
	}
	if (use.INLA.prev) prev.file.prefix <- paste0(prev.file.prefix, "INLA_")
    if (exclude.eldest.prev) prev.file.prefix <- paste0(prev.file.prefix, "no_elderly_")
    prev.mean.files <- paste0(prev.file.prefix, regions, "ons_meanlogprev2.txt")
    prev.sd.files <- paste0(prev.file.prefix, regions, "ons_sdlogprev2.txt")
    prev.file.prefix <- paste0(data.dirs["prev"], "/date_prev_", prev.file.txt, "_INLA_")
    prev.dat.file <- paste0(prev.file.prefix, "ons_dat2.csv")
    if(!format.inputs) names(prev.mean.files) <- names(prev.sd.files) <- regions
} else {
    prev.mean.files <- NULL
    prev.sd.files <- NULL
}
if(vacc.flag){
    vac1.files <- file.path(data.dirs["vacc"], paste0(str.date.vacc, "_1stvaccinations_", regions, ".txt"))

    if(vac.n_doses == 3) {
       vac2.files <- file.path(data.dirs["vacc"], paste0(str.date.vacc,"_2ndvaccinations_", regions, ".txt"))
       vac3.files <- file.path(data.dirs["vacc"], paste0(str.date.vacc,"_3rdvaccinations_", regions, ".txt"))
    } else {
        vacn.files <- file.path(data.dirs["vacc"], paste0(str.date.vacc,"_nthvaccinations_", regions, ".txt")) 
    }
#print("got here")
    if(!format.inputs) {    
    if(vac.n_doses == 3) {                                                                                                                                                   
        vac2.files <- file.path(data.dirs["vacc"], paste0(str.date.vacc,"_2ndvaccinations_", regions, ".txt"))                                                                    
        vac3.files <- file.path(data.dirs["vacc"], paste0(str.date.vacc,"_3rdvaccinations_", regions, ".txt"))

        load(build.data.filepath(file.path("RTM_format", region.type, "vaccination"), paste0(region.type, "vacc", str.date.vacc, ".RData")))                                 
        names(vac1.files) <- names(vac2.files)<- names(vac3.files) <- regions  

    } else {
        vacn.files <- file.path(data.dirs["vacc"], paste0(str.date.vacc,"_nthvaccinations_", regions, ".txt"))
        
        load(build.data.filepath(file.path("RTM_format", region.type, "vaccination"), paste0(region.type, "vacc", str.date.vacc, ".RData")))                            
        names(vac1.files) <- names(vacn.files) <- regions  
    }
#    print("got here 2")
    }
} else vac1.files <- vacn.files <- vac2.files <- vac3.files <- NULL
if(adm.flag){
    if(!format.inputs) {
        adm_csv_fname <- ifelse(admissions_only.flag, paste0("admissions_data_admissions_only", ifelse(cutoff_hosps_early & !deaths.flag & !hosp.flag, paste0("_drophosp_", gsub("-", "",toString(date_early_cutoff_hosps))), ""),".csv"),
                                                     paste0("admissions_data_all_hosp", ifelse(cutoff_hosps_early & !deaths.flag & !hosp.flag, paste0("_drophosp_", gsub("-", "",toString(date_early_cutoff_hosps))), ""),".csv"))
        adm.sam <- read_csv(file.path(data.dirs["adm"], adm_csv_fname))
        file.copy(file.path(data.dirs["adm"], adm_csv_fname), out.dir)
        file.rename(file.path(out.dir, adm_csv_fname), file.path(out.dir, "admissions_data.csv"))
        admsam.files <- paste0(data.dirs["adm"], "/", date.adm.str, "_", regions, "_", nA_adm, "ag_counts",ifelse(admissions_only.flag & data.desc == "admissions", "_adm_only", ""), ifelse(cutoff_hosps_early & !deaths.flag & !hosp.flag, paste0("_drophosp_", gsub("-", "",toString(date_early_cutoff_hosps))), ""),".txt")
    }
}


if(format.inputs){

    if(deaths.flag){
        if(data.desc == "reports") {
            source(file.path(proj.dir, "R/data/format_death_reports.R"))
        } else if (grepl("adjusted", data.desc)) {
            source(file.path(proj.dir, "R/data/format_adjusted_deaths.R"))
        } else if (running.England) {
            source(file.path(proj.dir, "R/data/format_deaths.R"))
        }
        if ("Scotland" %in% regions) {
            source(file.path(proj.dir, "R/data/format_Scottish_deaths.R"))
        }
        if ("Northern_Ireland" %in% regions) {
            source(file.path(proj.dir, "R/data/format_ni_deaths.R"))
        }
        if ("Wales" %in% regions) {
            source(file.path(proj.dir, "R/data/format_wales_deaths.R"))
        }
    }
    if(adm.flag) {
        source(file.path(proj.dir, "R/data/format_hosp_admissions.R"))
    }
    if(prev.flag){
        source(file.path(proj.dir, "R", "data", "format_prev.R"))
    }
    if(sero.flag){    ## Setup serology inputs
        source(file.path(proj.dir, "R/data/format_sero.R"))
    }
    if(vacc.flag){
        source(file.path(proj.dir, "R", "data", "format_vaccinations.R"))
    print('format_vacc.R sourced')
    }
    if(gp.flag){
        source(file.path(proj.dir, "R/data/format_linelist.R"))
    }
} 

print('got here')
## Set up the model specification.
source(file.path(proj.dir, "set_up.R"))
print('got here 2')
## Compile the code
if(compile.code) {
    system("make clean")
    system("make rtm_hpc2")
}

## Set up a requisite number of chains
startwd <- getwd()
setwd(out.dir)
if(exists("outpp"))
    rm(outpp)

cat(c(out.dir,'here'),'\n')
save.image("tmp.RData")
out.dir.orig <- out.dir

if(use.previous.run.for.start & length(previous.run.to.use) > 1){
    nchains <- length(previous.run.to.use)
    for(ichain in 2:nchains){
        new.dir <- paste0(out.dir, "_chain", ichain)
        cat(c(new.dir,'here_new'),'\n')
        if(dir.exists(new.dir))
            unlink(new.dir, recursive = TRUE)
        R.utils::copyDirectory(out.dir, new.dir)
        previous.loc <- previous.run.to.use[ichain]
        source(file.path(proj.dir, "import_pars.R"))
        source(file.path(proj.dir, "par_check.R"))
        knit(input = pars.template.loc, output = file.path(new.dir, "mod_pars.txt"))
        knit(input = inputs.template.loc, output = file.path(new.dir, "mod_inputs.txt"))
        out.dir <- new.dir
        setwd(out.dir)
        save.image("tmp.RData")
        out.dir <- out.dir.orig
        }
    }

setwd(out.dir)
if(run.code){
    system(file.path(proj.dir, "rtm_optim"), intern = TRUE)
	 system("chmod a-w coda* NNI* posterior* adaptive*")
}

## Post processing the results.
Rfile.loc <- file.path(file.loc, "R/output")
if(run.outputs){
    source(file.path(Rfile.loc, "tracePlots.R"))
	external = FALSE
	render(
		file.path(Rfile.loc, 'report-updated.Rmd'),
		html_document(pandoc_args = "--self-contained"),
		output_dir = out.dir,
		clean = FALSE, intermediates_dir = out.dir
	)
}

## Return back to initial directory
setwd(startwd)

