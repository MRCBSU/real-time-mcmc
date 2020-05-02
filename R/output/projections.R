## regions <- c("London", "North", "South.East", "South.West", "Midlands.and.East")
## region.index <- 4
## reg <- regions[region.index]
## wk <- 4

## negbin <- file.exists("coda_negbin_overdispersion")
## outputfile <- "./projection"
## datafile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_counts.txt")
## denomfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_denoms.txt")
## dataprojfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_counts.txt")
## denomprojfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_denoms.txt")

## time.horizon <- 39
## nproc <- 3

## nitr <- 10000

###### WHERE IS THE PROJECT ROUTE DIRECTORY
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
if(!exists("Rfile.loc")){
    Rfile.loc <- dirname(thisFile())
    proj.dir <- dirname(dirname(Rfile.loc))
    source(file.path(proj.dir, "set_up_inputs.R"))
}

###### WHERE IS THE R FILE DIRECTORY
target.dir <- out.dir
source(file.path(Rfile.loc, "input_extraction_fns.R"))

###### DIRECTORY CONTAINING MCMC OUTPUT
mcmc.file <- file.path(target.dir, "mcmc.RData")
if (!file.exists(mcmc.file)){
    source(file.path(Rfile.loc, "tracePlots.R"))
} else if(!exists("params")) load(mcmc.file)

## Last day of data
nt <- max(end.gp, end.hosp)

## start.date <- as.Date("02/10/2017", format = "%d/%m/%Y")
## dates <- start.date + 0:(time.horizon - 1)
## weeks <- format(dates, "%Y%V")

## population <- c(8787892, 15282212, 9026297, 5515953, 16655713)[region.index]
## data <- read.table(datafile, sep = "\t", row.names = 1)
## denoms <- read.table(denomfile, sep = "\t", row.names = 1) / population
## denoms <- denoms[, , drop = T]
## data <- data[, , drop = T]
## nt <- length(denoms)
## na <- 1
## nreg <- 5

## source("./generate.mu.for.pres.R")
iter.sum <- seq(from = i.saved / i.summary, to = i.saved, length.out = i.summary)

q.NNI <- list()
for(reg in regions){
    ## q.GP.ma <- apply(GP, 2, quantile, probs = c(0.025, 0.5, 0.975))
    q.NNI[[reg]] <- apply(NNI[[reg]], 2:3, sum)
    q.NNI[[reg]] <- apply(q.NNI[[reg]], 1, quantile, probs = c(0.025, 0.5, 0.975))

    pdf(file.path(target.dir, paste0("NNI_projections_", reg, ".pdf")))
    
    plot(dates.used, q.NNI[[reg]][2, ], type = "l", main = paste("Reconstructed (daily) Infections -", reg), ylab = "New Infections", xlab = "Day", ylim = c(0, max(q.NNI[[reg]])))
    lines(dates.used, q.NNI[[reg]][1, ], lty = 3)
    lines(dates.used, q.NNI[[reg]][3, ], lty = 3)
    abline(v = dates.used[nt], col = "red")
    dev.off()

    pdf(file.path(target.dir, paste0("log_NNI_projections_", reg, ".pdf")))
    
    plot(dates.used, log(q.NNI[[reg]][2, ]), type = "l", main = paste("Reconstructed (daily) Infections -", reg), ylab = "New Infections", xlab = "Day", ylim = c(0, log(max(q.NNI[[reg]]))))
    lines(dates.used, log(q.NNI[[reg]][1, ]), lty = 3)
    lines(dates.used, log(q.NNI[[reg]][3, ]), lty = 3)
    abline(v = dates.used[nt], col = "red")
    dev.off()

}


Rt.func <- function(vecS, matM){
    if(length(vecS) != nrow(matM)) stop("Dimension mismatch between vecS and matM")
    M.star <- sweep(matM, 2, vecS, `*`)
    max(abs(eigen(M.star)$value))
}

Rt <- q.Rt <- list()
M.star <- M <- M.mult <- list()
m <- params$contact_parameters[iter.sum, ]
R0 <- posterior.R0[iter.sum, , drop = F]
for(idir in 1:length(cm.bases)){
    M[[idir]] <- as.matrix(read_tsv(cm.bases[idir], col_names = FALSE))
    M.mult[[idir]] <- as.matrix(read_tsv(cm.mults[idir], col_names = FALSE)) + 1
    M.star[[idir]] <- array(apply(m, 1, function(mm) M[[idir]] * mm[M.mult[[idir]]]),
                            dim = c(nA, nA, nrow(m)))
}
m.levels <- cut(1:ndays, c(0, cm.breaks, Inf))
names(M) <- names(M.mult) <- names(M.star)
pop.total <- all.pop[1, ]
for(reg in regions){
    ireg <- which(regions %in% reg)
    M.temp <- matrix(M.star[[1]][, , 1, drop = F], dim(M.star[[1]])[1], dim(M.star[[1]])[2])
    R.star <- Rt.func(regions.total.population[ireg, ] / pop.total[ireg], M.temp)
    S <- apply(NNI[[reg]], c(1, 3), cumsum)  ## TxAxI array
    S <- -sweep(S, 2, regions.total.population[ireg, ], `-`) ## TxAxI
    R.prime <- sapply(1:ndays,
                      function(x) sapply(1:i.summary, function(i){
                          M.temp <- matrix(M.star[[m.levels[x]]][, , i, drop = F], nA, nA)
                          Rt.func(S[x, , i] / pop.total[ireg], M.temp)
                      }
                      )
                      ) ## IxT array
    Rt[[reg]] <- R0[, ireg] * R.prime / R.star ## I*T array
    q.Rt[[reg]] <- apply(Rt[[reg]], 2, quantile, probs = c(0.025, 0.5, 0.975))
}
## ############################################################################ ##

## ### Deaths ### ##
ifr <- params$prop_case_to_hosp[iter.sum, , drop = F]
source(file.path(Rfile.loc, "convolution.R"))
source(file.path(Rfile.loc, "gamma_fns.R"))
ifh <- (3*ifr) + 0.02
ifi <- 0.3 * ifh
## delay.to.hosp <- list(
##     incub.mean = 4,
##     incub.sd = 1.41,
##     disease.mean = 0,
##     disease.sd = 0,
##     report.mean = 7.21,
##     report.sd = 2.17)
## F.hosp <- discretised.delay.cdf(delay.to.hosp, steps.per.day = 1)
## delay.to.ICU <- list(
##                     incub.mean = 4,
##                     incub.sd = 1.41,
##                     disease.mean = 0,
##                     disease.sd = 0,
##                     report.mean = 9.21,
##                     report.sd = 3.59
##                     )
## F.icu <- discretised.delay.cdf(delay.to.ICU, steps.per.day = 1)
delay.to.death <- list(
    incub.mean = 4,
    incub.sd = 1.41,
    disease.mean = ddelay.mean,
    disease.sd = ddelay.sd,
    report.mean = rdelay.mean,
    report.sd = rdelay.sd
    )
F.death <- discretised.delay.cdf(delay.to.death, steps.per.day = 1)

## Load in some data for goodness-of-fit
dth.dat <- dth.dat.agg <- list()
if(all(file.exists(data.files))){
    for(r in 1:nr){
        dl <- end.hosp + reporting.delay
        dth.dat[[r]] <- read_tsv(data.files[r], col_names = FALSE)
        dth.dat[[r]] <- cbind(dth.dat[[r]][, 1], dth.dat[[r]][, -1])
        dth.dat.agg[[r]] <- cbind(dth.dat[[r]][, 1], apply(dth.dat[[r]][, -1], 1, sum))
        names(dth.dat)[r] <- names(dth.dat.agg)[r] <- regions[r]
    }
}

## Set up lists to store convolutions
H <- list()
ICU <- list()
D <- D.agg <- list()
rD <- list()
q.H <- list()
q.ICU <- list()
q.D <- list()
q.D.tot <- list()
q.rD <- list()
q.rD.tot <- list()
for(reg in regions){
    ireg <- which(regions %in% reg)
    ## Merge the two youngest age categories
    if(nA > 1)
        NNI[[reg]] <- apply(NNI[[reg]], 2:3, function(x) c(x[1]+x[2], x[-(1:2)]))
    ## H[[reg]] <- apply(NNI[[reg]], 2, function(x) x * t(ifh))
    ## H[[reg]] <- apply(H[[reg]], 1, conv, b = F.hosp)[1:(dim(H[[reg]])[2]), , drop = F]
    ## ICU[[reg]] <- apply(NNI[[reg]], 2, function(x) x * t(ifi))
    ## ICU[[reg]] <- apply(ICU[[reg]], 1, conv, b = F.icu)[1:(dim(ICU[[reg]])[2]), , drop = F]
    D[[reg]] <- apply(NNI[[reg]], c(1, 3), conv, b = F.death)[1:ndays, , , drop = F]
    D[[reg]] <- apply(D[[reg]], 1, function(x) x * t(ifr))
    D[[reg]] <- array(D[[reg]], dim = c(dim(t(ifr)), as.integer(ndays)))
    
    ## ICU[[reg]] <- apply(ICU[[reg]], c(1, 3), sum)
    ## D[[reg]] <- apply(D[[reg]], c(1, 3), sum)
    
    ## q.H[[reg]] <- apply(H[[reg]], 1, quantile, probs = c(0.025, 0.5, 0.975))
    ## q.ICU[[reg]] <- apply(ICU[[reg]], 1, quantile, probs = c(0.025, 0.5, 0.975))

    ## Get total deaths over age groups
    rD[[reg]] <- aperm(D[[reg]], c(2:3, 1)) ## I * T * A array
    q.D[[reg]] <- apply(D[[reg]], c(1, 3), quantile, probs = c(0.025, 0.5, 0.975))
    D.agg[[reg]] <- apply(D[[reg]], 2:3, sum)
    q.D.tot[[reg]] <- apply(D.agg[[reg]], 2, quantile, probs = c(0.025, 0.5, 0.975))

    ## pdf(file.path(target.dir, paste0("Hosp_projections_", reg, ".pdf")))
    ## plot(dates.used, q.H[[reg]][2, ], type = "l", main = paste("Projected (daily) Hospital Admissions - ", reg), ylab = "New Admissions", xlab = "Day", ylim = c(0, max(q.H[[reg]])))
    ## lines(dates.used, q.H[[reg]][1, ], lty = 3)
    ## lines(dates.used, q.H[[reg]][3, ], lty = 3)
    ## abline(v = dates.used[nt], col = "red")
    ## dev.off()

    ## pdf(file.path(target.dir, paste0("ICU_projections_", reg, ".pdf")))

    ## plot(dates.used, q.ICU[[reg]][2, ], type = "l", main = paste("Projected (daily) ICU Admissions - ", reg), ylab = "New Admissions", xlab = "Day", ylim = c(0, max(q.ICU[[reg]])))
    ## lines(dates.used, q.ICU[[reg]][1, ], lty = 3)
    ## lines(dates.used, q.ICU[[reg]][3, ], lty = 3)
    ## abline(v = dates.used[nt], col = "red")
    ## dev.off()

    ## Goodness-of-fit
    rD[[reg]] <- rnbinom(length(rD[[reg]]),
                         size = rD[[reg]] / params$hosp_negbin_overdispersion[iter.sum],
                         mu = rD[[reg]])
    rD[[reg]] <- array(rD[[reg]], dim = c(i.summary, ndays, ifelse(nA > 1, nA - 1, 1)))
    q.rD[[reg]] <- apply(rD[[reg]], 2:3, quantile, probs = c(0.025, 0.5, 0.975))
    q.rD.tot[[reg]] <- apply(apply(rD[[reg]], 1:2, sum), 2, quantile, probs = c(0.025, 0.5, 0.975))

    pdf(file.path(target.dir, paste0("Deaths_projections_", reg, ".pdf")), width = 20, height = 7)

    plot(dates.used, q.D.tot[[reg]][2, ], type = "l", main = paste("Projected (daily) Deaths - ", reg), ylab = "Count", xlab = "Day", ylim = c(0, max(q.rD.tot[[reg]])))
    lines(dates.used, q.D.tot[[reg]][1, ], lty = 3)
    lines(dates.used, q.D.tot[[reg]][3, ], lty = 3)
    abline(v = dates.used[nt], col = "red")

    points(dth.dat.agg[[reg]][start.hosp:dl, 1],
           dth.dat.agg[[reg]][start.hosp:dl, 2],
           pch = 18,
           col = "blue")
    arrows(dth.dat.agg[[reg]][start.hosp:dl, 1] + 0.1,
           q.rD.tot[[reg]][1, start.hosp:dl],
           y1 = q.rD.tot[[reg]][3, start.hosp:dl],
           code = 3,
           angle = 90,
           length = 0.1)
    
    dev.off()

    for(ag in 2:(nA - 1)){
        pdf(file.path(target.dir, paste0("Deaths_projections_", reg, "_age_", ag, ".pdf")),
            width = 20,
            height = 7)
        
        plot(dates.used,
             q.D[[reg]][3, ag, ],
             type = "l",
             main = paste("Projected (daily) Deaths -", reg, "among", age.labs[ag + 1], "year-olds"),
             ylab = "Count",
             xlab = "Day",
             lty = 3,
             ylim = c(0, max(q.rD[[reg]][3, , ag])))
        lines(dates.used, q.D[[reg]][1, ag, ], lty = 3)
        lines(dates.used, q.D[[reg]][2, ag, ])
        abline(v = dates.used[nt], col = "red")
        points(dth.dat[[reg]][start.hosp:dl, 1],
               dth.dat[[reg]][start.hosp:dl, ag + 2],
               pch = 18,
               col = "blue")
        arrows(dth.dat[[reg]][start.hosp:dl, 1] + 0.1,
               q.rD[[reg]][1, start.hosp:dl, ag],
               y1 = q.rD[[reg]][3, start.hosp:dl, ag],
               code = 3,
               angle = 90,
               length = 0.1)
        
        ## plot(dates.used, q.
        dev.off()
    }
}

## ## Get weekly tables of quantiles
require(dplyr)
require(tidyr)

get.tab <- function(x, dates){
    nni <- t(x)
    colnames(nni) <- dates
    nni <- as.data.frame(nni)
    nni <- nni %>%
        mutate(iter = 1:nrow(nni)) %>%
        pivot_longer(cols = -(nrow(x)+1), names_to = "Date", values_to = "Count")
    nni$Date <- as.Date(as.integer(nni$Date), origin = "1970-01-01")
    nni <- nni %>%
        mutate(Week = format(Date, "%Y-%W"))
    nni <- aggregate(nni$Count, by = list(iter = nni$iter, week = nni$Week), sum)
    aggregate(nni$x, by = list(Week = nni$week), quantile, probs = c(0.025, 0.5, 0.975))
}

nni <- hosp <- icu <- deaths <- NULL
for(reg in regions){
    NNI.sum <- apply(NNI[[reg]], 2:3, sum)
    nni <- rbind(nni, get.tab(NNI.sum, dates.used) %>%
                      mutate(Region = reg))
    ## hosp <- rbind(hosp, get.tab(H[[reg]], dates.used) %>%
    ##                     mutate(Region = reg))
    ## icu <- rbind(icu, get.tab(ICU[[reg]], dates.used) %>%
    ##                   mutate(Region = reg))
    deaths <- rbind(deaths, get.tab(t(D.agg[[reg]]), dates.used) %>%
                            mutate(Region = reg))
}

## weeks <- gl(n = ncol(GP.ma) / 7, k = 7, length = ncol(GP.ma))

## GP.ma.weeks <- by(t(GP.ma), weeks, function(x) apply(x, 2, sum))

## data.proj <- read.table(dataprojfile, sep = "\t", row.names = 1)
## denoms.proj <- read.table(denomprojfile, sep = "\t", row.names = 1) / population
## denoms.proj <- denoms.proj[, , drop = T]
## data.proj <- data.proj[, , drop = T]
## nt.proj <- length(denoms.proj)
## if(nt == nt.proj) nt <- nt - 1

## if(negbin){
##     mu.GP.ma <- t(t(GP[, 1:nt.proj]) * c(denoms[1:nt], denoms.proj[(nt + 1):nt.proj]))
##     r <- mu.GP.ma / mcmc$eta[, region.index]
##     X.GP <- array(rnbinom(length(r), size = r, mu = mu.GP.ma), dim = dim(r))
##     X.GP[is.na(X.GP)] <- 0
##     ## X.GP <- sapply(GP.ma, function(x){
##     ##     r <- x / inputs$eta
##     ##     rnbinom(length(x), size = r, mu = x)
##     ## })
## } else {
##     X.GP <- rpois(length(GP[, 1:nt]), GP[, 1:nt])
## }

## ## q.GP.ma.weeks <- sapply(GP.ma.weeks, quantile, probs = c(0.025, 0.5, 0.975))
## q.X.GP <- apply(X.GP, 2, quantile, probs = c(0.025, 0.5, 0.975))

## qGP <- round(quantile(mcmc$pgp, probs = c(0.025, 0.5, 0.975)), digits = 2)
## aR <- apply(NNI, 1, sum)
## qaR <- quantile(aR, probs = c(0.025, 0.5, 0.975))
## qaR <- round((qaR / population) * 100, digits = 4)

## pdf(paste0("gof_wash_", reg.pdf, "_wk", wk, ".pdf"))

## plot(q.GP.ma[2, 1:nt.proj] * c(denoms[1:nt], denoms.proj[(nt + 1):nt.proj]), type = "l", ylim = c(0, max(q.GP.ma[, 1:nt.proj] * denoms.proj, data.proj)), main = paste0("Goodness-of-fit to GP data - ", reg, "\n pGP = ", qGP[2], " (", qGP[1], ", ", qGP[3], "); AR% = ", qaR[2], " (", qaR[1], ", ", qaR[3], ")"), ylab = "Consultations", xlab = "Day",xlim = c(0, nt.proj))
## lines(q.GP.ma[1, 1:nt.proj] * c(denoms[1:nt], denoms.proj[(nt + 1):nt.proj]), lty = 3)
## lines(q.GP.ma[3, 1:nt.proj] * c(denoms[1:nt], denoms.proj[(nt + 1):nt.proj]), lty = 3)


## for(i in 1:nt){
##     arrows(i, y0 = q.X.GP[1, i], y1 = q.X.GP[3, i], angle = 90, code = 3, col = "yellowgreen", length = 0.1)
##     }
## for(i in (nt+1):nt.proj){
##     plotcol <- ifelse((nt+1) == nt.proj, "yellowgreen", "darkgreen")
##     arrows(i, y0 = q.X.GP[1, i], y1 = q.X.GP[3, i], angle = 90, code = 3, col = plotcol, length = 0.1)
##     }
## points(data, pch = 19, col = "violetred")
## plotcol <- ifelse((nt+1) == nt.proj, "violetred", "blue")
## points((nt+1):nt.proj, data.proj[(nt+1):nt.proj], pch = 19, col = plotcol)

## dev.off()

save(## q.ICU,
     q.D,
     q.NNI,
    q.Rt,
    q.rD,
     dates.used, file = file.path(target.dir, "plotted_summaries.RData"))

save(nni,
     ## icu,
     ## hosp,
     deaths, file = file.path(target.dir, "table_summaries.RData"))
