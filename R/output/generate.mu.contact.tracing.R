require(coda)
cmdArgs <- commandArgs(trailingOnly = FALSE)
scenario.number <- cmdArgs[1]
scenario.number <- 1
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else if (.Platform$GUI == "RStudio" || Sys.getenv("RSTUDIO") == "1") {
                # We're in RStudio
                return(rstudioapi::getSourceEditorContext()$path)
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

###### WHERE IS THE PROJECT ROUTE DIRECTORY
if(!exists("Rfile.loc")) Rfile.loc <- dirname(thisFile())
proj.dir <- dirname(dirname(Rfile.loc))
out.dir <- file.path(proj.dir, "model_runs", "20200501", "spim_20200501")

## out.dir <- "./"
## if(!exists("mcmc")) print("No mcmc object loaded\n")
## if(!exists("outputfile")) outputfile <- paste(out.dir, "../Xin/paper/mcmc/MCMC_andre_negbin/projections", sep = "")
## if(!exists("time.horizon")) time.horizon <- 329
## if(!exists("nproc")) nproc <- 1

## Function to read in MCMC arrays
read.mcmc <- function(bin.file, niter, npar){
    xtmp <- readBin(bin.file, "double", niter * npar)
    array(xtmp, dim = c(npar, niter))
}

## How many iterations, days, age groups, regions, etc...
n.iters <- 1025
n.proc <- 3
regions <- c("East_of_England", "London", "Midlands", "North_East_and_Yorkshire", "North_West", "South_East", "South_West")
n.reg <- length(regions)
time.horizon <- 440
n.age <- 8

## Read in the binary coda files.
mcmc.files <- list(hosp.disp = file(file.path(out.dir, "coda_hosp_negbin_overdispersion"), "rb"),
                   m = file(file.path(out.dir, "coda_contact_parameters"), "rb"),
                   egr = file(file.path(out.dir, "coda_exponential_growth_rate"), "rb"),
                   ip = file(file.path(out.dir, "coda_infectious_period"), "rb"),
                   nu = file(file.path(out.dir, "coda_log_p_lambda_0"), "rb"),
                   p.hosp = file(file.path(out.dir, "coda_prop_case_to_hosp"), "rb")
                   )
mcmc <- list(hosp.disp = read.mcmc(mcmc.files$hosp.disp, n.iters, 1),
             m = read.mcmc(mcmc.files$m, n.iters, 2 * n.reg),
             egr = read.mcmc(mcmc.files$egr, n.iters, n.reg),
             ip = read.mcmc(mcmc.files$ip, n.iters, 1),
             nu = read.mcmc(mcmc.files$nu, n.iters, n.reg),
             p.hosp = read.mcmc(mcmc.files$p.hosp, n.iters, n.age - 1),
             theta = 2/3
             )

source(file.path(Rfile.loc, "gamma_fns.R"))
source(file.path(Rfile.loc, "seir_reporting_contact_tracing_model_functions.R"))
source(file.path(Rfile.loc, "convolution.R"))
### control parameters
intervals.per.day <- 2

## population size
population <- c(74103, 318183, 804260, 704025, 1634429, 1697206, 683583, 577399, 122401, 493480, 1123981, 1028009, 3063113, 2017884, 575433, 483780, 118454, 505611, 1284743, 1308343, 2631847, 2708355, 1085138, 895188, 92626, 396683, 1014827, 1056588, 2115517, 2253914, 904953, 731817, 79977, 340962, 851539, 845215, 1786666, 1820128, 710782, 577678, 94823, 412569, 1085466, 1022927, 2175624, 2338669, 921243, 801040, 55450, 241405, 633169, 644122, 1304392, 1496240, 666261, 564958)
population <- matrix(population, n.reg, n.age, byrow = TRUE)
rownames(population) <- regions

## Set up vectors for key calculated parameters
R0 <- I0.tot <- vector("numeric", n.reg)
## Some other parameters we'll use
## Probability of being infectious before contact is traced
p <- 1 - (0.2^0.25) ## the probability of a geometric random variable giving a reading of < 4 time units
N <- rgeom(1000000, p)
X1 <- rgamma(1000000, shape = 2, rate = 1/2) ## Sample some latent periods.
idx <- which(N > 2*X1)
SEIR.param.list <- list(susceptibility = 1, ### A VECTOR OF INITIAL SUSCEPTIBILITIES BY AGE
                        R0.peak.day = 355,
                        relative.susceptibility.I1.to.I2 = 1,
                        latent.period = 4,
                        test.sensitivity = 1, ### THIS IS A DEFAULT VALUE, IF THERE IS AN ESTIMATED TEST.SENSITIVITY, THIS DOES NOT NEED TO BE CHANGED, ONLY IF A DIFFERENT FIXED VALUE IS USED
                        p.E = 0.8
                        )
## Times from infection to death
delay.to.death <- list(
    incub.mean = 4,
    incub.sd = 1.41,
    disease.mean = 15,
    disease.sd = 12.1,
    report.mean = 0,
    report.sd = 0
)

M.breakpoint <- c(36, 43, 50, 57, 64, 85, 106, 136, 181, 199, 251, 258, 307, 321, 363, 370, 410, 426)
M.intervals <- vector("list", length(M.breakpoint) + 1)
for(i in 1:length(M.intervals)){
  start.range <- ifelse(i == 1, 1, (M.breakpoint[i - 1] * intervals.per.day) + 1)
  end.range <- intervals.per.day * ifelse(i == length(M.intervals), time.horizon, M.breakpoint[i])
  M.intervals[[i]] <- start.range:end.range
}
M.files <- file.path(proj.dir, "contact_mats", c("england_8ag_contact.txt", "england_8ag_contact_ldwk1_20200501.txt", "england_8ag_contact_ldwk2_20200501.txt", "england_8ag_contact_ldwk3_20200501.txt", "england_8ag_contact_ldwk4_20200501.txt", "england_8ag_contact_ldwk5_20200501.txt", paste0("england_8ag_contact_scenario", scenario.number, "phase", c("1", "2", "3", "4a", "4b"), "_20200501.txt")))
M.param.files <- file.path(proj.dir, "contact_mats", paste0("ag8_mult", 0:5, ".txt"))[c(1, rep(2, 5), 3:5, rep(6, 10))]
M.files <- c(M.files, rep(M.files[10:11], 4))
get.contact.matrix <- function(Mfile, scale = FALSE)
    matrix(scan(Mfile, ), n.age, n.age) / (ifelse(scale, 1e7, 1))
M <- lapply(M.files, get.contact.matrix, scale = TRUE)
M.param <- lapply(M.param.files, get.contact.matrix)

mixing.model <- lapply(regions, function(reg){
    NGM <- population[reg, ] * M[[1]]
    evec <- eigen(NGM)$vectors[, 1]
    list(intervals = M.intervals,
         base.contact.matrices = M,
         scalants.param = M.param,
         initial.age.distribution = evec / sum(evec))
})
names(mixing.model) <- regions

## ###### Are we accounting for day of the week effects?
## read.design.matrix <- function(np, target.dir, design.file)
##     {
##         if(!is.null(design.file)){
##             X.mat <- scan(paste(target.dir, design.file, sep = ""), )
##             X.mat <- t(matrix(X.mat, np, length(X.mat) / np))
##         } else {
##             X.mat <- diag(np)
##         }
##         X.mat
##     } #### HERE!!!

## day.of.week.effect <- TRUE
## if(day.of.week.effect) dow.design.file <- "../Model4_GP_and_Viro_regional/d_o_w_proj_design_file.txt"
## dow.r.breaks <- FALSE

## dow.mat <- read.design.matrix(ncol(mcmc$dow),
##                               "./",
##                               dow.design.file)

## ### Day of week effect model
## if(!exists("region.index")) region.index <- 1
## if(day.of.week.effect)
## {
##   dow.effects <- list(design.matrix = dow.mat
##                       )
##   if(dow.r.breaks)
##     {
##       design.indices <- ((region.index - 1) * nrow(dow.mat)) + 1:nrow(dow.mat)
##       dow.effects$design.matrix <- dow.effects$design.matrix[design.indices, ]
##     }
## } else dow.effects <- NULL

## ## Proportion infected who go to ICU
## pICU <- rnorm(nitr, log(0.0001195139), sd = 0.51)
## pICU <- exp(pICU)

## ## Background consultation model
## bg.matrix <- read.design.matrix(ncol(mcmc$bg) / length(regions), "../Background_HHH/", "bgp_design_file.txt")
## library(Matrix)
## bg.matrix <- bdiag(bg.matrix, bg.matrix, bg.matrix, bg.matrix, bg.matrix)
## bg.matrix <- matrix(as.numeric(bg.matrix), nrow(bg.matrix), ncol(bg.matrix))
## bg.t.breaks <- seq(7, time.horizon, by = 7)
## bg.a.breaks <- NULL
## bg.T <- length(bg.t.breaks) + 1
## bg.A <- length(bg.a.breaks) + 1
## par.rows <- bg.T * bg.A
## bg.design.rows <- ((region.index - 1) * par.rows) + (1:par.rows)
## par.matrix <- matrix(1:par.rows, bg.A, bg.T)
## adf <- diff(c(0, bg.a.breaks, na))
## background.gp.model <- list(
##   regional.heterogeneity = TRUE,
##   breakpoints = bg.t.breaks,
##   design.matrix = bg.matrix[bg.design.rows, ],
##   parameterisation.matrix = par.matrix[rep(1:nrow(par.matrix), adf), , drop = F]
##   )


get.params <- function(mcmc, iter, ir){
    list(hosp.disp = mcmc$hosp.disp[, iter],
         mm = mcmc$m[(2*ir)-(1:0), iter],
         egr = mcmc$egr[ir, iter],
         aip = 2 + mcmc$ip[, iter],
         nu = mcmc$nu[ir, iter],
         p.hosp = mcmc$p.hosp[, iter],
         theta = mcmc$theta,
         pgp = 1)
}

reg <- regions[1]
ir <- which(regions %in% reg)

phi <- get.params(mcmc, 1, 1)

SEEIIR.model.function(phi = phi,
                      SEIR.param.list = SEIR.param.list,
                      mixing.model = mixing.model[[reg]],
                      population.size = population[reg, ],
                      country.population = sum(population),
                      num.days = time.horizon,
                      chasing.threshold = scenario.number * 15)

## ## if(!file.exists(paste("Plots_", outputfile, sep = "")))
## ##   system(paste("mkdir Plots_", outputfile, sep = ""))
## require(multicore)
## proc.block <- nr %/% nproc
## proc.remainder <- nr %% nproc
## proc.blocks <- c(rep(proc.block + 1, proc.remainder), rep(proc.block, nproc - proc.remainder))
## cumsum.blocks <- c(0, cumsum(proc.blocks))
## proc.args <- list(block.size = proc.block, remainder = proc.remainder, blocks = proc.blocks, cum.blocks = cumsum.blocks)
## procs <- list()
## for(proc in 1:nproc){
##     procs[[proc]] <- parallel({

##         GP <- NNI <- array(0, dim = c(proc.args$blocks[proc], time.horizon))
##         GP.gof <- array(0, dim = c(proc.args$blocks[proc], nrow(dow.mat)))
##         tdm <- dow.effects
##         tbm <- background.gp.model
##         for(inti in 1:proc.args$blocks[proc])
##         {
##             i <- proc.args$cum.blocks[proc] + inti

##             if(!is.null(tdm))
##                 tdm$params <- mcmc$dow[i, ]

##             mm <- mcmc$m[i, ]
##             mixing.model$M <- mixing.model$base.contact.matrices
            
##             for(j in 1:length(mixing.model$base.contact.matrices))
##                 mixing.model$M[[j]] <- mixing.model$base.contact.matrices[[j]] * matrix(mm[mixing.model$scalants.param[[j]]], nrow(mixing.model$scalants.param[[j]]), ncol(mixing.model$scalants.param[[j]]))
        
##             mixing.model$scaled.matrices <- scaled.mixing.matrix(mixing.model, population)
##             attach(mcmc)            
##             R0[i] <- egr[i, region.index] * aip[i] * (((egr[i, region.index] * SEIR.param.list$latent.period / 2) + 1)^2) / (1 - (1 / (((egr[i, region.index] * aip[i] / 2) + 1)^2)))
##             I0.tot[i] <- exp(nu[i, region.index]) * aip[i] * population / (pgp[i] * R0[i])
##             tmp <- SEEIIR.model.function(egr = egr[i, region.index],
##                                          aip = aip[i],
##                                          I0.tot = I0.tot[i],
##                                          SEIR.param.list = SEIR.param.list,
##                                          mixing.model = mixing.model,
##                                          num.days = time.horizon,
##                                          population.size = population)
            
##             ## NNI.ma <- SEEIIR.mass.action.model.function(egr = egr[i],
##             ##                                             aip = aip[i],
##             ##                                             I0.tot = I0.tot[i],
##             ##                                             SEIR.param.list = SEIR.param.list,
##             ##                                             mixing.model = mixing.model,
##             ##                                             num.days = time.horizon,
##             ##                                             population.size = population)$per.timestep

##             if(exists("pgp.model")){
##                 q <- c(pgp[i], pgp2[i])
##             } else {
##                 q <- pgp[i]
##                 pgp.model <- NULL
##             }
            
##             tbm$basic.rates <- bg[i, ]

##             NNI[inti, ] <- tmp$per.day
##             GP[inti, ] <- reporting.model(tmp$per.timestep,
##                                        q = q,
##                                        proportion.model = pgp.model,
##                                        background.model = tbm,
##                                        proportion.symptomatics = theta[i],
##                                        reporting.param.list = delay.to.GP,
##                                        num.days = time.horizon,
##                                        dow.effects = tdm,
##                                        population.size = population)$total
            
##             ## GP.gof[inti, ] <- reporting.model(NNI,
##             ##                               q = q,
##             ##                               proportion.symptomatics = theta[i],
##             ##                               proportion.model = pgp.model,
##             ##                               background.model = tbm,
##             ##                               dow.effects = tdm,
##             ##                               reporting.param.list = delay.to.GP,
##             ##                               num.days = nrow(dow.mat),
##             ##                               population.size = population)$total

##             ## ICU[inti, ] <- reporting.model(NNI,
##             ##                                   q = pICU[i],
##             ##                                   proportion.symptomatics = 1,
##             ##                                   proportion.model = pgp.model,
##             ##                                   reporting.param.list = delay.to.ICU,
##             ##                                   num.days = time.horizon,
##             ##                                   population.size = population)$total

##             ## pdf(paste("Plots_", outputfile, "/outputs_", i, ".pdf", sep = ""))
            
##             ## plot(GP[i, ], type = "l", ylim = c(0, max(GP[i, ], GP.ma[i, ])))
##             ## lines(GP.ma[i, ], lty = 2, col = "green")
##             detach(mcmc)
##             ## dev.off()
##             if(is.null(pgp.model)) rm(pgp.model)
##         }
##         list(GP = GP, NNI = NNI, R0 = R0)
##     })
## }
## res <- collect(procs)

## GP <- NNI <- R0 <- NULL
## for(i in 1:nproc)
##     {
##         GP <- rbind(GP,  res[[i]]$GP)
##         ## GP.gof <- rbind(GP.gof, res[[i]]$GP.gof)
##         ## ICU <- rbind(ICU, res[[i]]$ICU)
##         NNI <- rbind(NNI, res[[i]]$NNI)
##         R0 <- c(R0, res[[i]]$R0)
##     }

## write.table(GP, paste(outputfile, ".txt", sep = ""), row.names = F, col.names = F)
## ## write.table(GP.ma, paste(outputfile, "MASSACTION.txt", sep = ""), row.names = F, col.names = F)

## R0 <- R0[R0 != 0]
## qR0 <- quantile(R0, probs = c(0.025, 0.5, 0.975))
## qR0 <- round(qR0, digits = 4)

## save.image(paste(outputfile, ".RData", sep = ""))

