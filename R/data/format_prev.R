#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################

## Where to find the data as downloaded from the ONS?
## --------------------------------------------------
## Either already specified in the variable prev.loc
## or found as the latest entry into a default
## directory. If prev.loc is NULL, look for a
## command-line argument.
if(!exists("prev.loc")){ ## Set to default format for the filename
    input.loc <- "data/raw/prevalence"
    ## List the possible files in the directory
    prev.loc <- file.info(file.path(input.loc,
                                    list.files(path = input.loc,
                                               pattern=glob2rx(paste("202", tolower(region.type), "csv", sep = "*")))
                                    )
                          )
    ## Pick the most recently added
    input.loc <- rownames(prev.loc)[which.max(prev.loc$ctime)]
} else {
    if(is.null(prev.loc)){
        input.loc <- commandArgs(trailingOnly = TRUE)[1]
    } else {
        if(startsWith(prev.loc, "/")) input.loc <- prev.loc
        else input.loc <- build.data.filepath(subdir = "raw", "prevalence", prev.loc)
    }
}

## What is the date of publication of these data? If not specified, try to extract from filename
if(date.prev != lubridate::ymd(strapplyc(input.loc, "[0-9/]{8,}", simplify = TRUE))) stop("Specified date.prev does not match the most recent prevalence data file.")
## Substitute this date into the output file names
prev.mean.files <- gsub("date_prev", date.prev, prev.mean.files, fixed = TRUE)
prev.sd.files <- gsub("date_prev", date.prev, prev.sd.files, fixed = TRUE)

## Define an age-grouping
if(!exists("age.agg")){
    age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
    age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+")
    }
nA <- length(age.labs)

## Map our names for columns (LHS) to data column names (RHS)
possible.col.names <- list(
    age = c("age", "age_rtm", "ageg_rtm"),
    region = c("nhsregion", "onsregion"),
    sample_date = "date",
    day = "study_day",
    lmean = "lmean",
    lsd = "lsd"
    )
input.col.names <- suppressMessages(names(read_csv(input.loc, n_max=0)))
is.valid.col.name <- function(name) {name %in% input.col.names}
first.valid.col.name <- function(names) {first.where.true(names, is.valid.col.name)}
col.names <- lapply(possible.col.names, first.valid.col.name)
invalid.col.names <- sapply(col.names, is.null)
if (any(invalid.col.names)) {
	names.invalid.cols <- paste0(names(possible.col.names)[invalid.col.names], collapse = ", ")
	stop(paste("No valid column name for:", names.invalid.cols))
}

## Given a row in the prev data file, return its region, formatted with no spaces
get.region <- function(x) {
	if (region.type == "NHS") north.east.name <- "North_East_and_Yorkshire"
	if (region.type == "ONS") north.east.name <- "North_East"
    x %>% mutate(region = str_replace_all(region, " ", "_"),
				 region = recode(region,
					North_East = north.east.name,
					East_England = "East_of_England",
					North_East_England = north.east.name,
					North_West_England = "North_West",
					South_East_England = "South_East",
					South_West_England = "South_West",
					Yorkshire = "Yorkshire_and_The_Humber",
				))
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
if(!exists("prev.mean.files")){
    prev.mean.files <- build.data.filepath("RTM_format/prevalence",
                                           "prev",
                                           date.prev,
                                           "_",
                                           regions,
                                           "_ons_meanlogprev.txt")
    prev.sd.files <- build.data.filepath("RTM_format/prevalence",
                                         "prev",
                                         date.prev,
                                         "_",
                                         regions,
                                         "_ons_sdlogprev.txt")
}

## Which columns are we interested in?
prev.col.args <- list()
prev.col.args[[col.names[["sample_date"]]]] <- col_character()
prev.col.args[[col.names[["age"]]]] <- col_factor(ordered = TRUE)
prev.col.args[[col.names[["region"]]]] <- col_character()
prev.col.args[[col.names[["lmean"]]]] <- col_double()
prev.col.args[[col.names[["lsd"]]]] <- col_double()
prev.col.args[[col.names[["day"]]]] <- col_integer()
prev.cols <- do.call(cols, prev.col.args)

## Reading in the data ##
print(paste("Reading from", input.loc))
## strPos <- c("+", "Positive", "positive")
prev.dat <- read_csv(input.loc,
                     col_types = prev.cols) %>%
    rename(!!!col.names) %>%
	mutate(sample_date = fuzzy_date_parse(sample_date))
levels(prev.dat$age) <- age.labs[-1]

should.include.rows <- function(x) {
	include.days.mask <- x$day %in% prev.lik.days
	excluded.ages <- c("1-4")
	if (exclude.eldest.prev) excluded.ages <- c(excluded.ages, "65-74", "75+")
	include.ages.mask <- !(x$age %in% excluded.ages)
	return(include.days.mask & include.ages.mask)
}

## Filter to only those on or before the last day used in the likelihood.
prev.dat <- prev.dat %>%
    filter(sample_date <= (start.date + max(prev.lik.days) - 1)) %>%
    mutate(day = sample_date - start.date + 1) %>%
    ## Set equal to zero all those entries that are not going to be used in the likelihood.
    mutate(include = should.include.rows(.),
           lmean = ifelse(include, lmean, 0),
           lsd = ifelse(include, lsd, 0)) %>%
    select(-include) %>%
    get.region() %>%
    filter(region %in% regions)

# Check correct number of days
stopifnot(
  (
    prev.dat %>%
      filter(lsd > 0) %>%
      distinct(day, region) %>%
	  nrow()
  )
  == length(prev.lik.days) * nr
)

## Pad data with some zeros to take it back to day 1.
if(min(prev.dat$day) > 1){
    add.dates <- lubridate::as_date(start.date:(min(prev.dat$sample_date) - 1)) %>%
        as_tibble() %>%
        rename(sample_date = value)
    for(ag in age.labs)
        add.dates <- mutate(add.dates, !!ag := 0)
}

## Get into format for use by the rtm
names(prev.mean.files) <- names(prev.sd.files) <- regions
for(reg in regions){
    region.mean <- bind_rows(add.dates,
                             pivot_wider(prev.dat %>%
                               filter(region == reg),
                               id_cols = sample_date,
                               names_from = age,
                               values_from = lmean) %>%
                             mutate(!!age.labs[1] := 0)
                             )
    region.sd <- bind_rows(add.dates,
                           pivot_wider(prev.dat %>%
                               filter(region == reg),
                               id_cols = sample_date,
                               names_from = age,
                               values_from = lsd) %>%
                             mutate(!!age.labs[1] := 0)
                           )

    print(paste("Writing to",
                prev.mean.files[reg],
                "."))

    region.mean %>%
        write_tsv(prev.mean.files[reg], col_names = FALSE)

    print(paste("Writing to",
                prev.sd.files[reg],
                "."))

    region.sd %>%
        write_tsv(prev.sd.files[reg], col_names = FALSE)
	stopifnot(max(region.mean[-1]) > 0)
	stopifnot(max(region.sd[-1]) > 0)
    
}
write_csv(prev.dat, prev.dat.file)
