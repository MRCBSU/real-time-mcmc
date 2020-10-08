source('run.R')
setwd(out.dir)
if(file.exists(file.path(out.dir, "coda_lfx"))) stop("Output already exists")
system('sbatch ../../../submit_cu_hpc_with_postprocessing')
