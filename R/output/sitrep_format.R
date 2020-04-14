require(dplyr)
require(tidyr)
require(tidyverse)
require(lubridate)

## regions <- c("London", "North", "South.East", "South.West", "Midlands.and.East")
## region.index <- 4
## reg <- regions[region.index]
## wk <- 4
rtm.dir <- proj.dir <- "./"
R.dir <- paste0(rtm.dir, "R/output/")
res.dir <- paste0(rtm.dir, "model_runs/20200412/", c("L_OL_report_variable_relax/", "L_OL_report_variable_tight/"))
run.date <- lubridate::as_date("20200412")

str.scenario <- c("forecast - relaxed prior", "forecast - tight prior")
death_type <- "death_inc_line"
## Dates of forecast period
first.date <- run.date - 21
end.date <- run.date + 21
day.zero <- as.Date("16/02/2020", format = "%d/%m/%Y")

## negbin <- file.exists("coda_negbin_overdispersion")
## outputfile <- "./projection"
## datafile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_counts.txt")
## denomfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_denoms.txt")
## dataprojfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_counts.txt")
## denomprojfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_denoms.txt")

## time.horizon <- 39
## nproc <- 3

## nitr <- 10000

source(paste0(R.dir, "/convolution.R"))
source(paste0(R.dir, "/gamma_fns.R"))

out.all <- NULL

for(dir in res.dir){
    
    load(paste0(dir, "mcmc.RData"))
    idir <- which(res.dir %in% dir)
    
    ## Last day of data
    
    ## population <- c(8787892, 15282212, 9026297, 5515953, 16655713)[region.index]
    ## data <- read.table(datafile, sep = "\t", row.names = 1)
    ## denoms <- read.table(denomfile, sep = "\t", row.names = 1) / population
    ## denoms <- denoms[, , drop = T]
    ## data <- data[, , drop = T]
    ## nt <- length(denoms)
    ## na <- 1
    ## nreg <- 5

    ifr <- params$prop_case_to_hosp[seq(10, nrow(params$hosp_negbin_overdispersion), length.out = 1000), , drop = F]
    ifh <- (3*ifr) + 0.02
    ifi <- 0.3*ifh
    delay.to.hosp <- list(
        incub.mean = 4,
        incub.sd = 1.41,
        disease.mean = 0,
        disease.sd = 0,
        report.mean = 7.21,
        report.sd = 2.17)
    F.hosp <- discretised.delay.cdf(delay.to.hosp, steps.per.day = 1)
    delay.to.ICU <- list(
        incub.mean = 4,
        incub.sd = 1.41,
        disease.mean = 0,
        disease.sd = 0,
        report.mean = 9.21,
        report.sd = 3.59
    )
    F.icu <- discretised.delay.cdf(delay.to.ICU, steps.per.day = 1)
    delay.to.death <- list(
        incub.mean = 4,
        incub.sd = 1.41,
        disease.mean = 0,
        disease.sd = 0,
        report.mean = 17.8,
        report.sd = 8.9
    )
    F.death <- discretised.delay.cdf(delay.to.death, steps.per.day = 1)
    
    ## Function to format output tables
    get.tab <- function(x, dates, reg, scenario, type, run.date, start.date, end.date, weekflag = FALSE, qprobs = c(0.5, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)){
        nni <- t(x)
        colnames(nni) <- dates
        nni <- as.data.frame(nni)
        nni <- nni %>%
            mutate(iter = 1:nrow(nni)) %>%
            pivot_longer(cols = -(nrow(x)+1), names_to = "Date", values_to = "Value")
        nni$Date <- as.Date(as.integer(nni$Date), origin = "1970-01-01")
        if(weekflag){
            nni <- nni %>%
                mutate(Week = format(Date, "%Y-%W"))
            nni <- aggregate(nni$Value, by = list(iter = nni$iter, week = nni$Week), sum)
            nni <- aggregate(nni$x, by = list(Week = nni$week), quantile, probs = qprobs)
        } else nni <- aggregate(nni$Value, by = list(Date = nni$Date), quantile, probs = qprobs)
        ## qstrings <- paste0("x.", as.character(100*qprobs), "%")
        nni %>%
            mutate(Group = "PHE",
                   CreationDate = run.date,
                   Scenario = scenario,
                   Geography = reg,
                   ValueType = type) %>%
            rename(Value = x) %>%
            select(Group,
                   Scenario,
                   CreationDate,
                   Date,
                   Geography,
                   ValueType,
                   everything()) %>%
            filter(Date > start.date,
                   Date <= end.date)
    }
    
    
    names(regions.total.population) <- regions
    new.regions <- list(England = regions)
    ## new.regions <- NULL
    regions.total.population <- c(regions.total.population, sum(regions.total.population[new.regions[[1]]]))
    names(regions.total.population) <- c(regions, names(new.regions))
    ## regions <- c(
    ## 	"East_of_England",
    ## 	"London",
    ## 	"Midlands",
    ## 	"North_East_and_Yorkshire",
    ## 	"North_West",
    ## 	"South_East",
    ## 	"South_West"
    ## )
    
    NNI.out <- list() ## What we're going to put into the report
    AR <- list()
    H <- list()
    ICU <- list()
    D <- list()
    ar <- nni <- hosp <- icu <- deaths <- NULL
    for(reg in regions){
        
        ## Initially working with no age stratification.
        NNI.out[[reg]] <- apply(NNI[[reg]][1,,], 2, cumsum)
                                        #AR[[reg]] <- apply(NNI.out[[reg]], 2, cumsum) / regions.total.population[reg]
                                        #AR[[reg]] <- AR[[reg]][run.date - day.zero, , drop = F]
        
        H[[reg]] <- apply(NNI[[reg]], 2, function(x) x * t(ifh))
        H[[reg]] <- apply(H[[reg]], 1, conv, b = F.hosp)[1:(dim(H[[reg]])[2]), , drop = F]
        ICU[[reg]] <- apply(NNI[[reg]], 2, function(x) x * t(ifi))
        ICU[[reg]] <- apply(ICU[[reg]], 1, conv, b = F.icu)[1:(dim(ICU[[reg]])[2]), , drop = F]
        D[[reg]] <- apply(NNI[[reg]], c(1, 3), conv, b = F.death)[1:(dim(NNI[[reg]])[2]), , , drop = F]
        D[[reg]] <- apply(D[[reg]], 1, function(x) x * t(ifr))
        D[[reg]] <- t(D[[reg]]) ## NEEDS SOME CONSIDERATION OF AGE
    }
    for(r in 1:length(new.regions)){
        reg <- names(new.regions)[r]
        NNI[[reg]] <- do.call('+', NNI[new.regions[[r]]])
        NNI.out[[reg]] <- do.call('+', NNI.out[new.regions[[r]]])
        ## AR[[reg]] <- apply(NNI.out[[reg]], 2, cumsum) / regions.total.population[reg]
        ## AR[[reg]] <- AR[[reg]][run.date - day.zero, , drop = F]
        H[[reg]] <- do.call('+', H[new.regions[[r]]])
        ICU[[reg]] <- do.call('+', ICU[new.regions[[r]]])
        D[[reg]] <- do.call('+', D[new.regions[[r]]])
        ## H[[reg]] <- apply(NNI.out[[reg]], 1, function(x) x * t(ifh))
        ## H[[reg]] <- apply(H[[reg]], 1, conv, b = F.hosp)[1:(dim(H[[reg]])[2]), , drop = F]
        ## ICU[[reg]] <- apply(NNI.out[[reg]], 1, function(x) x * t(ifi))
        ## ICU[[reg]] <- apply(ICU[[reg]], 1, conv, b = F.icu)[1:(dim(ICU[[reg]])[2]), , drop = F]
        ## D[[reg]] <- apply(NNI.out[[reg]], 2, conv, b = F.death)[1:(dim(NNI.out[[reg]])[1]), , drop = F]
        ## D[[reg]] <- apply(D[[reg]], 1, function(x) x * t(ifr))
        ## D[[reg]] <- t(D[[reg]])
    }
    for(reg in names(NNI.out)){
        inc <- rbind(nni,
                     get.tab(NNI[[reg]][1, , ],
                             dates.used,
                             reg,
                             str.scenario[idir],
                             type = "infections",
                             run.date = run.date,
                             start.date = first.date,
                             end.date = end.date)
                     )
        nni <- rbind(nni,
                     get.tab(NNI.out[[reg]],
                             dates.used,
                             reg,
                             str.scenario[idir],
                             type = "infections_cum",
                             run.date = run.date,
                             start.date = first.date,
                             end.date = end.date)
                     )
        deaths <- rbind(deaths,
                        get.tab(D[[reg]],
                                dates.used,
                                reg,
                                str.scenario[idir],
                                type = death_type,
                                run.date = run.date,
                                start.date = first.date,
                                end.date = end.date)
                        )
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
       
    out <- rbind(inc, nni, deaths) %>%
        filter(Geography == "England") %>%
        mutate(`1st centile` = Value[, 2],
               `5th centile` = Value[, 3],
               `25th centile` = Value[, 4],
               `50th centile` = Value[, 5],
               `75th centile` = Value[, 6],
               `95th centile` = Value[, 7],
               `99th centile` = Value[, 8],
               `Creation Day` = day(CreationDate),
               `Creation Month` = month(CreationDate),
               `Creation Year` = year(CreationDate),
               `Day of Value` = day(Date),
               `Month of Value` = month(Date),
               `Year of Value` = year(Date),
               Geography = str_replace_all(Geography, "_", " "),
               Model = str.scenario[idir])
    
    out$Value <- out$Value[, 1]
    
    out.all <- rbind(out.all, out %>% select(-c(CreationDate, Date)))
    
}

out.dir <- "sitep_outputs/"## , format(run.date, format = "%Y%m%d"))
if(!file.exists(out.dir))
    system(paste("mkdir", out.dir))
write.csv(out.all,
          paste0(out.dir, "sitrep", format(run.date, format = "%Y%m%d"), ".csv"),
          row.names = FALSE)
