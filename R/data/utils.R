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
