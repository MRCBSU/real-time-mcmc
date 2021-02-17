library(readr)
library(tidyr)
library(dplyr)

source("utils.R")
dir.data <- ("../../data/population")
## Load in the population data in long format
pop <- read_csv(build.data.filepath(subdir = "", "pop_long_nhsrlo_2018_syoa.csv"))
## Load in a lookup table to get the region names, rather than codes.
lookup <- read_csv(build.data.filepath(subdir = "", "nhs_region_lookup.csv"))
lookup <- unique(lookup[, c(3, 4)])
## If we're aggregating by age, get the aggregation being used or set to default.
if(!exists("age.agg")){
    age.agg <- c(0, 1, 5, 15, 25, 45, 65, 91)
    age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65+")
    }

## Put data into wide format - so it resembles the data sent by the ONS.
pop <- pop %>%
    filter(Sex == 4) %>%
    pivot_wider(id_cols = NHSERApr19CD, names_from = Age, values_from = Population)
pop <- pop %>%
    mutate("All ages" = apply(pop[, -1], 1, sum)) %>%
    inner_join(.,lookup) %>%
    rename(Code = NHSERApr19CD,
           Name = NHSERApr19NM) %>%
    mutate(Geography = "NHS Region") %>%
    select(Code,
           Name,
           Geography,
           "All ages",
           everything())

# Add the population of UK constituent nations
pop <- read_csv(build.data.filepath(subdir = "", "popn2018_all.csv")) %>%
	rename(Geography=Geography1) %>%
	filter(Geography=="Country") %>%
	bind_rows(pop) %>%
	assertr::assert(assertr::is_uniq, Name) %>%
       assertr::assert(function(x) !is.na(x), everything())

## Put back into long format with age aggregation.
## Long format...
regions <- pop %>%
    pivot_longer(-(1:4), names_to = "Age", values_to = "popn")
regions$Age <- as.integer(regions$Age)
regions$Grp <- cut(regions$Age, age.agg, age.labs,right=FALSE,ordered_result=T)
regions.short <- aggregate(regions$popn, by = list(Region=regions$Name,Age=regions$Grp), FUN = sum, drop=TRUE)
regions.short$Region <- gsub(" ", "_", regions.short$Region)


## Define some aggregate regions that may get used...
nhs.regions <- list("London" = "London",
                    "Outside_London" = c("North East and Yorkshire", "North West", "Midlands", "East of England", "South East", "South West"),
                    "UNITED_KINGDOM" = "UNITED KINGDOM",
                    "ENGLAND" = unlist(lookup[, 2]),
					"East_of_England" =	"East of England",
					"London" = "London",
					"Midlands" = "Midlands",
					"North_East_and_Yorkshire" = "North East and Yorkshire",
					"North_West" = "North West",
					"South_East" = "South East",
					"South_West" = "South West",
					"Scotland" = "SCOTLAND",
					"Wales" = "WALES",
					"Northern_Ireland" = "NORTHERN IRELAND"
                    )


save(pop, nhs.regions, file = build.data.filepath(subdir = "", "pop_nhs.RData"))
