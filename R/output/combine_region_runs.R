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
proj.dir <- dirname(dirname(file.loc))
source(file.path(proj.dir, "set_up_inputs.R"))

all.out.dirs <- file.path(proj.dir, "model_runs", subdir.name, all.regions)
names(all.out.dirs) <- all.regions

out.dir <- file.path(proj.dir, "model_runs", subdir.name, "_OVERALL_")
out.dir.correct <- out.dir
flg.createfile <- !file.exists(out.dir)
if(flg.createfile) system(paste("mkdir", out.dir))

all.q.NNI.cum <- list()
all.q.D.cum <- list()
all.q.NNI <- list()
all.q.D  <- list()
all.NNI <- list()
for (region in all.regions) {
	load(file.path(all.out.dirs[[region]], "occupancy_results.RData"))
	load(file.path(all.out.dirs[[region]], "plotted_summaries.RData"))
	load(file.path(all.out.dirs[[region]], "mcmc.RData"))
	all.q.NNI.cum[[region]] <- q.NNI.cum[[region]]
	all.q.D.cum[[region]] <- q.D.cum[[region]]
	all.q.NNI[[region]] <- q.NNI[[region]]
	all.q.D[[region]]  <- q.D[[region]]
	all.NNI[[region]] <- NNI[[region]]
	file.copy(file.path(all.out.dirs[[region]], "codas.pdf"),
				file.path(out.dir, paste0("codas_", region, ".pdf")))
}
q.NNI.cum <- all.q.NNI.cum
q.D.cum <- all.q.D.cum
q.NNI <- all.q.NNI
q.D  <- all.q.D
q.NNI  <- all.q.NNI
NNI <- all.NNI
out.dir <- out.dir.correct 
hosp.data <- build.data.filepath("RTM_format", "deaths", date.of.runs, "_", all.regions, ".txt")

source(file.path(file.loc, "build_report.R"))
