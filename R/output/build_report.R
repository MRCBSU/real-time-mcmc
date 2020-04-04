library(ggplot2)
library(lubridate)
library(tidyverse)
#knitr::opts_chunk$set(echo = FALSE)

start_date = ymd("2020-02-17")

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

day.number <- function(date) {
	return(time_length(date - start_date, "days") + 1)
}

file.loc <- dirname(thisFile())
proj.dir <- dirname(dirname(file.loc))

out.file <- function(...) {
	return(file.path(out.dir, paste0(..., collapse="")))
}
if (!exists("out.dir")) source(file.path(proj.dir, "set_up_inputs.R"))

if (!exists("q.NNI.cum")) {
	if (!file.exists(out.file("occupancy_results.RData"))) {
		source(file.path(file.loc, "icu_occupancy.R"))
	}
	load(out.file("mcmc.RData"))
	load(out.file("plotted_summaries.RData"))
	load(out.file("occupancy_results.RData"))
}

add.numerical.summary.for.day <- function(x, day) {
	day.num <- day.number(day)
	extract_quantile_value <- function(quantile) {
		return(sapply(x$quantile_list, function(y) {prettyNum(round(y[quantile, day.num]), big.mark=",")}))
	}
	args <- list(x)
	args[[format(day)]] <- paste0(extract_quantile_value("50%"), " (",
					 extract_quantile_value("2.5%"), "--",
					 extract_quantile_value("97.5%"), ")"
			)
	return(do.call(mutate, 	args))
}

numerical.summary <- tibble()
for (reg in names(q.NNI.cum)) {
	names <- paste0(
		c("Cumulative infections", "Cumulative deaths"),
		" (", reg, ")"
	)
	reg.summary <- tribble(
			~output_name,	~quantile_list,
			names[1],		q.NNI.cum[[reg]],
			names[2],		q.D.cum[[reg]],
		) %>%
		add.numerical.summary.for.day(today()) %>%
		add.numerical.summary.for.day(today() + 7) %>%
		add.numerical.summary.for.day(today() + 14)
	numerical.summary <- rbind(numerical.summary, reg.summary)
}


transform.q <- function(x, str = "Mid"){
	rownames(x) <- c("lower", "median", "upper")
    df <- as.data.frame(cbind(dates.used, t(x)))
	as_tibble(t(x)) %>% mutate(date = lubridate::as_date(dates.used))
}

plot.q <- function(q, ylab) {
	ggplot(data = q, aes(x = date, y = median, ymin = lower, ymax = upper)) +
		geom_line() +
		geom_ribbon(alpha=0.25) +
		geom_vline(xintercept=today(), linetype = "dashed", color = "red") +
		xlab("Date") +
		ylab(ylab) +
		theme_minimal()
}

names(hosp.data) <- names(NNI)


rmarkdown::render(file.path(file.loc, "report.Rmd"), "html_document", output_dir = out.dir)
