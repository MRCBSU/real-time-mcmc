library(ggplot2)
library(lubridate)
library(tidyverse)
#knitr::opts_chunk$set(echo = FALSE)

start_date = ymd("2020-02-17")
reg = "East_England"

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

out.file <- function(filename) {
	return(file.path(out.dir, filename))
}
source(file.path(proj.dir, "set_up_inputs.R"))

if (!file.exists(out.file("occupancy_results.RData"))) {
	source(file.path(file.loc, "icu_occupancy.R"))
}
load(out.file("plotted_summaries.RData"))
load(out.file("occupancy_results.RData"))

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

numerical.summary <- tribble(
	~output_name,				~quantile_list,
	"Cumulative infections",	q.NNI.cum,
	"Cumulative deaths",		q.D.cum,
	"Current ICU occupancy",	q.occupancy,
	) %>%
	add.numerical.summary.for.day(today()) %>%
	add.numerical.summary.for.day(today() + 7) %>%
	add.numerical.summary.for.day(today() + 14)


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

D.data <- read_tsv(hosp.data,  col_names = c("date", "incidence")) %>%
	mutate(date = ymd(date), cum = cumsum(incidence)) %>%
	filter(date >= ymd("2020-03-10"))
ll.data <- read_tsv(gp.data,  col_names = c("date", "incidence")) %>%
	mutate(date = ymd(date), cum = cumsum(incidence)) %>%
	filter(date >= ymd("2020-03-10"))
ll.data.plot <- ggplot(data = ll.data, aes(x = date, y = incidence)) +
	geom_col() +
	xlab("Date of lab report") +
	ylab("Confirmed cases")
ggsave(filename = out.file("confirmed_cases_data.png"), plot = ll.data.plot, width = 7, height = 15, unit = "cm")
D.data.plot <- ggplot(data = D.data, aes(x = date, y = incidence)) +
	geom_col() +
	xlab("Date of death") +
	ylab("Deaths reported")
ggsave(filename = out.file("deaths_data.png"), plot = D.data.plot, width = 7, height = 15, unit = "cm")


df.NNI <- transform.q(q.NNI$East)
plot.NNI <- plot.q(df.NNI, "New infections")
ggsave(filename = out.file("NNI.pdf"), plot = plot.NNI)

df.NNI.cum <- transform.q(q.NNI.cum)
plot.NNI.cum <- plot.q(df.NNI.cum, "Cumulative infections")
ggsave(filename = out.file("NNI_cum.pdf"), plot = plot.NNI.cum)

df.D <- transform.q(q.D$East_England) %>% left_join(D.data)
plot.D <- plot.q(df.D, "Incidence of deaths") + geom_point(aes(date, incidence))
ggsave(filename = out.file("deaths.pdf"), plot = plot.D)

df.D.cum <- transform.q(q.D.cum) %>% left_join(D.data)
plot.D.cum <- plot.q(df.D.cum, "Cumulative deaths") + geom_point(aes(date, cum))
ggsave(filename = out.file("D_cum.pdf"), plot = plot.D.cum)

df.icu.occ <- transform.q(q.occupancy)
plot.icu.occ <- plot.q(df.icu.occ, "ICU occupancy")
ggsave(filename = out.file("icu_cum.pdf"), plot = plot.icu.occ)

df.icu.cum <- transform.q(q.ICU$East_England)
plot.icu.cum <- plot.q(df.icu.cum, "ICU admissions")
ggsave(filename = out.file("icu_admission.pdf"), plot = plot.icu.cum)


rmarkdown::render(file.path(file.loc, "report.Rmd"), "html_document", output_dir = out.dir)
