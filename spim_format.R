require(dplyr)
require(tidyr)
require(tidyverse)
require(lubridate)

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
file.loc <- dirname(thisFile())
proj.dir <- file.loc
source(file.path(proj.dir, "config.R"))

## regions <- c("London", "North", "South.East", "South.West", "Midlands.and.East")
## region.index <- 4
## reg <- regions[region.index]
## wk <- 4
R.dir <- file.path(proj.dir, "R", "output")
res.dir <- combined.dir
run.date <- lubridate::ymd(date.data)

str.scenario <- "forecast"
first.date <- run.date - 30
end.date <- run.date + 21
## negbin <- file.exists("coda_negbin_overdispersion")
## outputfile <- "./projection"
## datafile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_counts.txt")
## denomfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_denoms.txt")
## dataprojfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_counts.txt")
## denomprojfile <- paste0("/project/pandemic_flu/Data/GP_In_Hours/8Mar18/", reg, "_denoms.txt")

## time.horizon <- 39
## nproc <- 3

## nitr <- 10000

load(file.path(res.dir, "mcmc.RData"))

## Last day of data
#nt <- 42
## Days of forecast

day.zero <- as.Date("16/02/2020", format = "%d/%m/%Y")
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

source(paste0(R.dir, "/convolution.R"))
source(paste0(R.dir, "/gamma_fns.R"))
## source("./generate.mu.for.pres.R")
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
get.tab <- function(x, dates, reg, scenario, type, run.date, start.date, end.date, weekflag = FALSE, qprobs = c(0.5, 0.01, 0.05, 0.25, 0.75, 0.95, 0.99)){
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
               Scenario = "Forecast",
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


#names(regions.total.population) <- regions
new.regions <- list(England = regions)
new.regions <- NULL
#regions.total.population <- c(regions.total.population, sum(regions.total.population[new.regions[[1]]]))
#names(regions.total.population) <- c(regions, names(new.regions))
regions <- c(
	"East_of_England",
	"London",
	"Midlands",
	"North_East_and_Yorkshire",
	"North_West",
	"South_East",
	"South_West",
	"Scotland"
)

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
for(reg in names(NNI.out)){
    nni <- rbind(nni,
                 get.tab(NNI.out[[reg]],
                         dates.used,
                         reg,
                         str.scenario,
                         type = "infections_cum",
                         run.date = run.date,
                         start.date = first.date,
                         end.date = end.date)
                 )
    deaths <- rbind(deaths,
                    get.tab(D[[reg]],
                            dates.used,
                            reg,
                            str.scenario,
                            type = "death_inc_line",
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

dir.string <- paste0("date_", run.date)
if(!file.exists(dir.string))
    system(paste("mkdir", dir.string))
## save(nni, icu, deaths, file = paste(dir.string, "table_summaries.RData", sep = "/"))
## write.csv(nni, paste0(dir.string, "/nni_", str.scenario, ".csv"))
## write.csv(icu, paste0(dir.string, "/icu_", str.scenario, ".csv"))
## write.csv(deaths, paste0(dir.string, "/deaths_", str.scenario, ".csv"))

out <- rbind(nni, deaths)
out <- out %>%
    mutate(`1st centile` = Value[, 2],
           `5th centile` = Value[, 3],
           `25th centile` = Value[, 4],
           `75th centile` = Value[, 5],
           `95th centile` = Value[, 6],
           `99th centile` = Value[, 7],
           `50th centile` = Value[, 7],
		   `Creation Day` = day(CreationDate),
		   `Creation Month` = month(CreationDate),
		   `Creation Year` = year(CreationDate),
		   `Day of Value` = day(Date),
		   `Month of Value` = month(Date),
		   `Year of Value` = year(Date),
		   Geography = str_replace_all(Geography, "_", " "))

out$Value <- out$Value[, 1]

out <- out %>% select(-c(CreationDate, Date))

write.csv(out, paste0(dir.string, "/PHEall_", str.scenario, "_v2.csv"), row.names = FALSE)

