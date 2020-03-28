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

icu.int <- round(ICU$East)
icu.out <- apply(icu.int, 2, row.leave.time)
icu.net <- icu.int - icu.out
icu.occupancy <- apply(icu.net, 2, cumsum)

pdf(file.path(target.dir, "ICU_occupancy.pdf"))
q.occupancy <- apply(icu.occupancy, 1, quantile, probs = c(0.025, 0.5, 0.975))
plot(dates.used, q.occupancy[2, ], type = "l", main = "ICU bed occupancy", ylab = "Beds occupied", xlab = "Day", ylim = c(0, max(q.occupancy)))
lines(dates.used, q.occupancy[1, ], lty = 3)
lines(dates.used, q.occupancy[3, ], lty = 3)
abline(v = dates.used[nt], col = "red")
dev.off()

pdf(file.path(target.dir, "NNI_cum.pdf"))
NNI.cum <- apply(NNI$East[1,,], 2, cumsum)
q.NNI.cum <- apply(NNI.cum, 1, quantile, probs = c(0.025, 0.5, 0.975))
plot(dates.used, q.NNI.cum[2, ], type = "l", main = "Cumulative infection count", ylab = "Cumulative infections", xlab = "Day", ylim = c(0, max(q.NNI.cum)))
lines(dates.used, q.NNI.cum[1, ], lty = 3)
lines(dates.used, q.NNI.cum[3, ], lty = 3)
abline(v = dates.used[nt], col = "red")
dev.off()

pdf(file.path(target.dir, "Deaths_cum.pdf"))
D.cum <- apply(D$East, 2, cumsum)
q.D.cum <- apply(D.cum, 1, quantile, probs = c(0.025, 0.5, 0.975))
plot(dates.used, q.D.cum[2, ], type = "l", main = "Cumulative deaths count", ylab = "Cumulative deaths", xlab = "Day", ylim = c(0, max(q.D.cum)))
lines(dates.used, q.D.cum[1, ], lty = 3)
lines(dates.used, q.D.cum[3, ], lty = 3)
abline(v = dates.used[nt], col = "red")
dev.off()

save(q.occupancy, q.NNI.cum, q.D.cum, file = file.path(target.dir, "occupancy_results.RData"))
