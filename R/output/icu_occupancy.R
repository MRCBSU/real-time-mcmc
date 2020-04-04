## Location of this script
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

## Where are various directories?
file.loc <- dirname(thisFile())

source(file.path(file.loc, 'projections.R'))

element.leave.time <- function(i, row) {
	result <- rep(0, length(row))
	for (j in seq(from=1, length=row[i])) {
		leave <- round(rgamma(1, 8, 1)) + i
		if (leave <= length(result)) {
			result[leave] <- result[leave] + 1
		}
	}
	return(result)
}

row.leave.time <- function(row) {
	return(Reduce('+',lapply(seq(from=1,to=length(row)), element.leave.time, row=row)))
}

if (FALSE) { # We need better method of estimating ICU occupancy and this is slow
pdf(file.path(target.dir, "ICU_occupancy.pdf"))
icu.int <- list()
icu.out <- list()
icu.net <- list()
icu.occupancy <- list()
q.occupancy <- list()
for (reg in names(NNI)) {
	icu.int[[reg]] <- round(ICU[[reg]])
	icu.out[[reg]] <- apply(icu.int[[reg]], 2, row.leave.time)
	icu.net[[reg]] <- icu.int[[reg]] - icu.out[[reg]]
	icu.occupancy[[reg]] <- apply(icu.net[[reg]], 2, cumsum)

	q.occupancy[[reg]] <- apply(icu.occupancy[[reg]], 1, quantile, probs = c(0.025, 0.5, 0.975))
	plot(dates.used, q.occupancy[[reg]][2, ], type = "l", main = "ICU bed occupancy", ylab = "Beds occupied", xlab = "Day", ylim = c(0, max(q.occupancy[[reg]])))
	lines(dates.used, q.occupancy[[reg]][1, ], lty = 3)
	lines(dates.used, q.occupancy[[reg]][3, ], lty = 3)
	abline(v = dates.used[nt], col = "red")
	dev.off()
}
}

NNI.cum = list()
q.NNI.cum = list()
pdf(file.path(target.dir, "NNI_cum.pdf"))
for (reg in names(NNI)) {
	NNI.cum[[reg]] <- apply(NNI[[reg]][1,,], 2, cumsum)
	q.NNI.cum[[reg]] <- apply(NNI.cum[[reg]], 1, quantile, probs = c(0.025, 0.5, 0.975))
	plot(dates.used, q.NNI.cum[[reg]][2, ], type = "l", main = "Cumulative infection count", ylab = "Cumulative infections", xlab = "Day", ylim = c(0, max(q.NNI.cum[[reg]])))
	lines(dates.used, q.NNI.cum[[reg]][1, ], lty = 3)
	lines(dates.used, q.NNI.cum[[reg]][3, ], lty = 3)
	abline(v = dates.used[nt], col = "red")
	dev.off()
}

D.cum = list()
q.D.cum = list()
pdf(file.path(target.dir, "Deaths_cum.pdf"))
for (reg in names(NNI)) {
	D.cum[[reg]] <- apply(D[[reg]], 2, cumsum)
	q.D.cum[[reg]] <- apply(D.cum[[reg]], 1, quantile, probs = c(0.025, 0.5, 0.975))
	plot(dates.used, q.D.cum[[reg]][2, ], type = "l", main = "Cumulative deaths count", ylab = "Cumulative deaths", xlab = "Day", ylim = c(0, max(q.D.cum[[reg]])))
	lines(dates.used, q.D.cum[[reg]][1, ], lty = 3)
	lines(dates.used, q.D.cum[[reg]][3, ], lty = 3)
	abline(v = dates.used[nt], col = "red")
	dev.off()
}

save(q.occupancy, q.NNI.cum, q.D.cum, file = file.path(target.dir, "occupancy_results.RData"))
