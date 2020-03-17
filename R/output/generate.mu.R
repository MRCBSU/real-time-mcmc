out.dir <- "./"
if(!exists("mcmc")) print("No mcmc object loaded\n")
if(!exists("outputfile")) outputfile <- paste(out.dir, "../Xin/paper/mcmc/MCMC_andre_negbin/projections", sep = "")
if(!exists("time.horizon")) time.horizon <- 329
if(!exists("nproc")) nproc <- 1

## OMITTED PARAMETERS?
nr <- nrow(mcmc$m)
if(length(grep("theta", names(mcmc))) == 0)
  mcmc$theta <- rep(0.15, nr)
if(length(grep("aip", names(mcmc))) == 0)
  mcmc$aip <- rep(3.45, nr)
if(!any(names(mcmc) == "m"))
    mcmc$m <- matrix(1, nr, 3)

homedir <- "~/"
#homedir <- "~/bsu_home/"
## homedir <- "/Volumes/Home_PaulB/"
proj.dir <- "/project/pandemic_flu/"
#proj.dir <- "~/bsu_pandemic/"
## proj.dir <- "/Volumes/Pandemic_flu/"

source(paste(homedir, "R/functions/statdist/gamma_fns.R", sep = ""))
source(paste(proj.dir, "Real_Time_Modelling/Rlib/seir_reporting_functions.R", sep = ""))

### control parameters
intervals.per.day <- 2

if(!exists("population")) population <- 55268067

R0 <- I0.tot <- vector("numeric", nr)
SEIR.param.list <- list(susceptibility = 0.625, ### A VECTOR OF INITIAL SUSCEPTIBILITIES BY AGE
                        R0.peak.day = 355,
                        relative.susceptibility.I1.to.I2 = 1,
                        latent.period = 2,
                        test.sensitivity = 1 ### THIS IS A DEFAULT VALUE, IF THERE IS AN ESTIMATED TEST.SENSITIVITY, THIS DOES NOT NEED TO BE CHANGED, ONLY IF A DIFFERENT FIXED VALUE IS USED
                        )
delay.to.GP <- list(
                    incub.mean = 1.62,
                    incub.sd = 1.791,
                    disease.mean = 1.968,
                    disease.sd = 1.199,
                    report.mean = 0.5,
                    report.sd = 0.5
                    )
delay.to.ICU <- list(
    incub.mean = 1.62,
    incub.sd = 1.791,
    disease.mean = 3.122093,
    disease.sd = 3.012615,
    report.mean = 0.5,
    report.sd = 0.5
)

M.breakpoint <- c(20, 27, 79, 93, 132, 139)
M.intervals <- vector("list", length(M.breakpoint) + 1)
for(i in 1:length(M.intervals)){
  start.range <- ifelse(i == 1, 1, (M.breakpoint[i - 1] * intervals.per.day) + 1)
  end.range <- intervals.per.day * ifelse(i == length(M.intervals), time.horizon, M.breakpoint[i])
  M.intervals[[i]] <- start.range:end.range
}
M <- vector("list", length(M.intervals))
M.param <- list(1, 3, 1, 3, 1, 3, 1)
for(i in 1:length(M)){
  M[[i]] <- matrix(1, 1, 1)
  M.param[[i]] <- matrix(M.param[[i]], 1, 1)
}


mixing.model <- list(initial.age.distribution = as.vector(1),
                     intervals = M.intervals,
                     base.contact.matrices = M,
                     scalants.param = M.param)
                     
###### Are we accounting for day of the week effects?
read.design.matrix <- function(np, target.dir, design.file)
    {
        if(!is.null(design.file)){
            X.mat <- scan(paste(target.dir, design.file, sep = ""), )
            X.mat <- t(matrix(X.mat, np, length(X.mat) / np))
        } else {
            X.mat <- diag(np)
        }
        X.mat
    } #### HERE!!!

day.of.week.effect <- TRUE
if(day.of.week.effect) dow.design.file <- "../Model4_GP_and_Viro_regional/d_o_w_proj_design_file.txt"
dow.r.breaks <- FALSE

dow.mat <- read.design.matrix(ncol(mcmc$dow),
                              "./",
                              dow.design.file)

### Day of week effect model
if(!exists("region.index")) region.index <- 1
if(day.of.week.effect)
{
  dow.effects <- list(design.matrix = dow.mat
                      )
  if(dow.r.breaks)
    {
      design.indices <- ((region.index - 1) * nrow(dow.mat)) + 1:nrow(dow.mat)
      dow.effects$design.matrix <- dow.effects$design.matrix[design.indices, ]
    }
} else dow.effects <- NULL

## Proportion infected who go to ICU
pICU <- rnorm(nitr, log(0.0001195139), sd = 0.51)
pICU <- exp(pICU)

## Background consultation model
bg.matrix <- read.design.matrix(ncol(mcmc$bg) / length(regions), "../Background_HHH/", "bgp_design_file.txt")
library(Matrix)
bg.matrix <- bdiag(bg.matrix, bg.matrix, bg.matrix, bg.matrix, bg.matrix)
bg.matrix <- matrix(as.numeric(bg.matrix), nrow(bg.matrix), ncol(bg.matrix))
bg.t.breaks <- seq(7, time.horizon, by = 7)
bg.a.breaks <- NULL
bg.T <- length(bg.t.breaks) + 1
bg.A <- length(bg.a.breaks) + 1
par.rows <- bg.T * bg.A
bg.design.rows <- ((region.index - 1) * par.rows) + (1:par.rows)
par.matrix <- matrix(1:par.rows, bg.A, bg.T)
adf <- diff(c(0, bg.a.breaks, na))
background.gp.model <- list(
  regional.heterogeneity = TRUE,
  breakpoints = bg.t.breaks,
  design.matrix = bg.matrix[bg.design.rows, ],
  parameterisation.matrix = par.matrix[rep(1:nrow(par.matrix), adf), , drop = F]
  )

## if(!file.exists(paste("Plots_", outputfile, sep = "")))
##   system(paste("mkdir Plots_", outputfile, sep = ""))
require(multicore)
proc.block <- nr %/% nproc
proc.remainder <- nr %% nproc
proc.blocks <- c(rep(proc.block + 1, proc.remainder), rep(proc.block, nproc - proc.remainder))
cumsum.blocks <- c(0, cumsum(proc.blocks))
proc.args <- list(block.size = proc.block, remainder = proc.remainder, blocks = proc.blocks, cum.blocks = cumsum.blocks)
procs <- list()
for(proc in 1:nproc){
    procs[[proc]] <- parallel({

        GP <- NNI <- array(0, dim = c(proc.args$blocks[proc], time.horizon))
        GP.gof <- array(0, dim = c(proc.args$blocks[proc], nrow(dow.mat)))
        tdm <- dow.effects
        tbm <- background.gp.model
        for(inti in 1:proc.args$blocks[proc])
        {
            i <- proc.args$cum.blocks[proc] + inti

            if(!is.null(tdm))
                tdm$params <- mcmc$dow[i, ]

            mm <- mcmc$m[i, ]
            mixing.model$M <- mixing.model$base.contact.matrices
            
            for(j in 1:length(mixing.model$base.contact.matrices))
                mixing.model$M[[j]] <- mixing.model$base.contact.matrices[[j]] * matrix(mm[mixing.model$scalants.param[[j]]], nrow(mixing.model$scalants.param[[j]]), ncol(mixing.model$scalants.param[[j]]))
        
            mixing.model$scaled.matrices <- scaled.mixing.matrix(mixing.model, population)
            attach(mcmc)            
            R0[i] <- egr[i, region.index] * aip[i] * (((egr[i, region.index] * SEIR.param.list$latent.period / 2) + 1)^2) / (1 - (1 / (((egr[i, region.index] * aip[i] / 2) + 1)^2)))
            I0.tot[i] <- exp(nu[i, region.index]) * aip[i] * population / (pgp[i] * R0[i])
            tmp <- SEEIIR.model.function(egr = egr[i, region.index],
                                         aip = aip[i],
                                         I0.tot = I0.tot[i],
                                         SEIR.param.list = SEIR.param.list,
                                         mixing.model = mixing.model,
                                         num.days = time.horizon,
                                         population.size = population)
            
            ## NNI.ma <- SEEIIR.mass.action.model.function(egr = egr[i],
            ##                                             aip = aip[i],
            ##                                             I0.tot = I0.tot[i],
            ##                                             SEIR.param.list = SEIR.param.list,
            ##                                             mixing.model = mixing.model,
            ##                                             num.days = time.horizon,
            ##                                             population.size = population)$per.timestep

            if(exists("pgp.model")){
                q <- c(pgp[i], pgp2[i])
            } else {
                q <- pgp[i]
                pgp.model <- NULL
            }
            
            tbm$basic.rates <- bg[i, ]

            NNI[inti, ] <- tmp$per.day
            GP[inti, ] <- reporting.model(tmp$per.timestep,
                                       q = q,
                                       proportion.model = pgp.model,
                                       background.model = tbm,
                                       proportion.symptomatics = theta[i],
                                       reporting.param.list = delay.to.GP,
                                       num.days = time.horizon,
                                       dow.effects = tdm,
                                       population.size = population)$total
            
            ## GP.gof[inti, ] <- reporting.model(NNI,
            ##                               q = q,
            ##                               proportion.symptomatics = theta[i],
            ##                               proportion.model = pgp.model,
            ##                               background.model = tbm,
            ##                               dow.effects = tdm,
            ##                               reporting.param.list = delay.to.GP,
            ##                               num.days = nrow(dow.mat),
            ##                               population.size = population)$total

            ## ICU[inti, ] <- reporting.model(NNI,
            ##                                   q = pICU[i],
            ##                                   proportion.symptomatics = 1,
            ##                                   proportion.model = pgp.model,
            ##                                   reporting.param.list = delay.to.ICU,
            ##                                   num.days = time.horizon,
            ##                                   population.size = population)$total

            ## pdf(paste("Plots_", outputfile, "/outputs_", i, ".pdf", sep = ""))
            
            ## plot(GP[i, ], type = "l", ylim = c(0, max(GP[i, ], GP.ma[i, ])))
            ## lines(GP.ma[i, ], lty = 2, col = "green")
            detach(mcmc)
            ## dev.off()
            if(is.null(pgp.model)) rm(pgp.model)
        }
        list(GP = GP, NNI = NNI, R0 = R0)
    })
}
res <- collect(procs)

GP <- NNI <- R0 <- NULL
for(i in 1:nproc)
    {
        GP <- rbind(GP,  res[[i]]$GP)
        ## GP.gof <- rbind(GP.gof, res[[i]]$GP.gof)
        ## ICU <- rbind(ICU, res[[i]]$ICU)
        NNI <- rbind(NNI, res[[i]]$NNI)
        R0 <- c(R0, res[[i]]$R0)
    }

write.table(GP, paste(outputfile, ".txt", sep = ""), row.names = F, col.names = F)
## write.table(GP.ma, paste(outputfile, "MASSACTION.txt", sep = ""), row.names = F, col.names = F)

R0 <- R0[R0 != 0]
qR0 <- quantile(R0, probs = c(0.025, 0.5, 0.975))
qR0 <- round(qR0, digits = 4)

save.image(paste(outputfile, ".RData", sep = ""))

