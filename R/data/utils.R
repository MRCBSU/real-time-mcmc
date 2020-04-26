## Given several parts of a filename, build the whole filepath
library(assertr)
build.data.filepath <- function(subdir, ...) {
	dir <- file.path(dir.data, subdir)
	filename <- paste0(...)
	return(gsub("//", "/", file.path(dir, filename)))
}


## Parse dates which might be in weird formats or even have times
fuzzy_date_parse <- function(date) {
	as_date(parse_date_time(date, c('dmyT', 'ymdT', 'dmyR', 'ymdR'), truncated=1))
}

swap.day.and.month <- function(date) {
	return(ymd(paste(year(date), day(date), month(date))))
}

## Functions governing the region mapping
mapping <- list(
	"East Midlands" = "Midlands",
	"East of England" = "East of England",
	"London" = "London",
	"North East" = "North East and Yorkshire",
	"North West" = "North West",
	"South East" = "South East",
	"South West" = "South West",
	"West Midlands" = "Midlands",
	"Yorkshire and Humber" = "North East and Yorkshire"
)
phe.to.nhs.region <- function(x) {
	args <- mapping
	args[[".x"]] <- x$phe_region
	result <- do.call(recode, args)
	result[x$utla_name == "Cumbria"] <- NA
	return(result)
}

nhs.region <- function(x) {
	# Below code will give the NHS region
	ifelse(
		!is.na(x$nhs_region),
		x$nhs_region,
		phe.to.nhs.region(x)
	) %>%
	str_replace_all(" ", "_")
}

inverse.map <- function(r, l){
    assign.NAs <- any(sapply(l, function(x) any(is.na(x))))
    if(is.na(r) & !assign.NAs) return("NA")
    bl <- sapply(l, function(x) r %in% x)
    ifelse(all(!bl), "NA", names(l)[bl])
}

## Vectorising inverse.map
map.to.region <- function(xcol)
    unlist(sapply(xcol, inverse.map, l = all.regions[regions]))

## Function for cleaning dates in data
within.range <- function(dates) {
	return(dates <= today() & dates >= ymd("2019-01-01"))
}

heuristically.swap.day.and.month <- function(dates) {
	swapped.dates <- suppressWarnings(swap.day.and.month(dates))
	should.swap <- (!is.na(swapped.dates) &
					!within.range(dates) &
					within.range(swapped.dates))
	dates[should.swap] <- swapped.dates[should.swap]
	return(dates)
}

## Region definitions
all.regions <- list( ## Region names that might be used and their constituent NHS regions.
    England = c("Midlands",
                "East_of_England",
                "London",
                "North_East_and_Yorkshire",
                "North_West",
                "South_East",
                "South_West",
                NA),
    East_of_England = "East_of_England",
    London = "London",
    Midlands = "Midlands",
    North_East_and_Yorkshire = "North_East_and_Yorkshire",
    North_West = "North_West",
    South_East = "South_East",
    South_West = "South_West"
)
all.regions$Outside_London = all.regions$England[-3]

