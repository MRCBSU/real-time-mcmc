library(ggplot2)
library(lubridate)
library(tidyverse)
#knitr::opts_chunk$set(echo = FALSE)

start_date <- ymd("2020-02-17")

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

calc.posterior.summary <- function(posterior, col=NULL) {
	if (ncol(posterior) > 1) posterior <- posterior[,col]
	quantiles <- quantile(posterior, c(0.025, 0.5, 0.975))
	return(tribble(
		~Median,		~`95% lower`,		~`95% upper`,
		quantiles[2],	quantiles[1],		quantiles[3]
	))
}

if(is.null(params$contact_parameters))
    params$contact_parameters <- t(array(contact.reduction, dim = dim(t(posterior.R0))))
    
if (is.null(names(posterior.R0))) {
    posterior.ifr <- params$prop_case_to_hosp
	posterior.summary <-
		calc.posterior.summary(posterior.R0) %>%
		bind_rows(calc.posterior.summary(params$contact_parameters)) %>%
		bind_cols(parameter = c("R0", "Lockdown effect"))
} else {
	R0.summary <- bind_rows(lapply(posterior.R0, calc.posterior.summary))
	col.names <- sapply(names(posterior.R0), function(x) {paste0("R0 (", str_replace_all(x, "_", " "), ")")})
	R0.summary$parameter <- col.names

	contact_param<- bind_rows(lapply(posterior.contact_param, calc.posterior.summary, col = 2))
	col.names <- sapply(names(posterior.contact_param), function(x) {paste0("Lockdown effect (", str_replace_all(x, "_", " "), ")")})
	contact_param$parameter <- col.names

	ifr <- bind_rows(lapply(posterior.ifr, calc.posterior.summary))
	col.names <- sapply(names(posterior.ifr), function(x) {paste0("IFR (", str_replace_all(x, "_", " "), ")")})
	ifr$parameter <- col.names

	num.today <- day.number(lubridate::ymd(date.data))
	Rt <- bind_rows(lapply(posterior.Rt, calc.posterior.summary, col = num.today))
	col.names <- sapply(names(posterior.ifr), function(x) {
							paste0("R on ", lubridate::today(), " (", str_replace_all(x, "_", " "), ")")
						})
	Rt$parameter <- col.names

	posterior.summary <- rbind(R0.summary, contact_param, ifr, Rt)
}

posterior.summary <- posterior.summary %>%
	select(parameter, everything()) %>%
	arrange(parameter)

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
predict.on <- lubridate::ymd(date.data)
for (reg in names(q.NNI.cum)) {
	names <- paste0(
		c("Cumulative infections", "Cumulative deaths", "Rt"),
		" (", reg, ")"
	)
	reg.summary <- tribble(
			~output_name,	~quantile_list,
			names[1],		q.NNI.cum[[reg]],
			names[2],		q.D.cum[[reg]],
		) %>%
		add.numerical.summary.for.day(predict.on) %>%
		add.numerical.summary.for.day(lubridate::ymd(20200413)) %>%
		add.numerical.summary.for.day(predict.on + 7) %>%
		add.numerical.summary.for.day(predict.on + 14)
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
