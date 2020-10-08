suppressMessages(library(lubridate))
suppressMessages(library(tidyverse))

source(file.path(proj.dir, "R/data/format_deaths.R"))
age_groupings <- c("0-44", "45-54", "55-64", "65-74", ">=75")
raw.deaths <- dth.dat %>%
    mutate(Age.Grp = cut(
		age,
		c(0, 45, 55, 65, 75, Inf),
		age_groupings,
		right = FALSE, ordered_result = T)) %>%
	count(Date, Region, Age.Grp)
#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

# Given a vector, return the first element where predicate is true
# Returns NULL if not true for any
first.where.true <- function(vec, predicate) {
	true.at <- which(predicate(vec))
	if (length(true.at) == 0) return(NULL)
	index.to.use <- min(true.at)
	return(vec[index.to.use])
}	

if(!exists("date.data"))
    date.data <- (today() - days(1)) %>% format("%Y%m%d")

## Where to find the data, if NULL use command line argument
possible.deaths.locations <- c(
	file.path(proj.dir, "data/raw/deaths/adjusted", paste0(ymd(date.data), ".csv")))
deaths.loc <- first.where.true(possible.deaths.locations, file.exists)
if (is.null(deaths.loc)) {
	stop(paste('No valid deaths data files, tried:', possible.deaths.locations))
}


####################################################################
## BELOW THIS LINE SHOULD NOT NEED EDITING
####################################################################

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
if(!exists("file.loc")){
    file.loc <- dirname(thisFile())
    proj.dir <- dirname(dirname(file.loc))
    dir.data <- file.path(proj.dir, "data")
	source(file.path(file.loc, "utils.R"))
	source(file.path(proj.dir, "config.R"))
}
if(!exists("data.files"))
    data.files <- build.data.filepath("RTM_format/deaths",
                                      "deaths",
                                      date.data,
                                      "_",
                                      regions,
                                      "_",
                                      nA,
                                      "ages.txt")

## Read the file and rename columns
latest.date <- ymd(date.data) ## - reporting.delay
earliest.date <- ymd("2020-02-17")

input.loc <- deaths.loc
print(paste("Reading from", input.loc))
dth.dat <- read_csv(
	input.loc,
	col_types = cols(
		region = col_character(),
		age_group = col_character(),
		onset_date = col_date(),
		posterior_mean = col_number(),
		posterior_median = col_integer(),
		posterior_low = col_number(),
		posterior_upp = col_number(),
		observed = col_integer()
	)
) %>%
	rename(Date = onset_date, n = posterior_median, Region = region) %>%
	mutate(
		Age.Grp = factor(
			age_group,
			levels = age_groupings,
			ordered = T
		),
		Region = recode(
			Region,
			east = "East_of_England",
			london = "London",
			midlands = "Midlands",
			northeast = "North_East_and_Yorkshire",
			northwest = "North_West",
			southeast = "South_East",
			southwest = "South_West"
		)
	) %>%
    filter(Date <= latest.date) %>%
    filter(Date >= earliest.date) %>%
    filter(Region %in% regions) %>%
	filter(Age.Grp %in% age_groupings, !is.na(age_group)) %>%
	full_join(
		raw.deaths,
		by = c("Region", "Date", "Age.Grp")
	) %>%
	right_join(		# Add missing rows
		expand_grid(
			Date = as_date(earliest.date:latest.date),
			Region = regions,
			Age.Grp = age_groupings
		),
		by = c("Date", "Region", "Age.Grp")
	) %>%
	mutate(n = ifelse(is.na(n.x), n.y, n.x),
		   Age.Grp = factor(Age.Grp, levels = age_groupings)) %>%
	replace_na(list(n = 0)) %>%
	arrange(Region, Date, Age.Grp) %>%
	select(Region, Date, Age.Grp, n)


## Write rtm.dat to data file
names(data.files) <- regions
for(reg in regions) {
    region.dat <- pivot_wider(dth.dat %>%
                              filter(Region == reg),
                              names_from = Age.Grp,
                              values_from = n) %>%
					select(-Region)
    tmpFile <- data.files[reg]
	dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)
    
    print(paste(
        "Writing to",
        tmpFile,
        "(", sum(region.dat[, -1]), "total deaths,", nrow(region.dat), "rows.)"
    ))
    
    region.dat %>%
        write_tsv(
            tmpFile,
            col_names = FALSE
        )
}

## Save the data as processed
save(dth.dat, rtm.dat, file = file.path(out.dir, "deaths_data.RData"))

## Save a quick plot of the data...
require(ggplot2)
rtm.dat %>%
    group_by(Date, Region) %>%
    summarise(count = sum(n)) %>%
    mutate(ignore = !(Date <= (latest.date - reporting.delay))) -> rtm.dat.plot

gp <- ggplot(rtm.dat.plot, aes(x = Date, y = count, color = Region)) +
    geom_line(aes(linetype = ignore)) +
    geom_point() +
    theme_minimal() +
    ggtitle(paste("Daily number of deaths by day of death (on", lubridate::as_date(date.data), ")")) +
    xlab("Date of death") +
    ylab("#Deaths") +
    theme(
        legend.position = "top",
        legend.justification = "left",
        )
plot.filename <- build.data.filepath("RTM_format/deaths", "deaths_plot", date.data, "_", reporting.delay, "d", ifelse(flg.cutoff, paste0("_cutoff", str.cutoff), ""), ".pdf")
if (!file.exists(dirname(plot.filename))) dir.create(dirname(plot.filename))
ggsave(plot.filename,
       gp + guides(linetype=FALSE),
       width = 1.5*8.5,
       height = 1.5*6)

## rtm.dat.plot %>%
##     group_by(Date, ignore) %>%
##     summarise(count = sum(count)) -> rtm.dat.Eng.plot
