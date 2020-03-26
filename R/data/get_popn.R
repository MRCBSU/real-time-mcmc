pop <- read.csv(build.data.filepath(subdir = "", "popn2018_all.csv"))
if(!exists("age.agg")){
    age.agg <- c(0, 1, 5, 15, 25, 45, 65, 91)
    age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65+")
    }

regions <- pop[pop$Geography1 == "Region", ]

diff(age.agg)

library(tidyr)

regions <- regions %>%
    pivot_longer(-(1:4), names_to = "Age", values_to = "popn")
regions$Age <- as.integer(gsub("X","",regions$Age,fixed=TRUE))
regions$Grp <- cut(regions$Age, age.agg, age.labs,right=FALSE,ordered_result=T)
regions$popn <- as.integer(gsub(",","",regions$popn, fixed = TRUE))
regions.short <- aggregate(regions$popn, by = list(Region=regions$Name,Age=regions$Grp), FUN = sum, drop=TRUE)
regions.short$Region <- regions.short$Region[, drop=TRUE]
## remove spaces from region names
levels(regions.short$Region) <- gsub(" ", "_", levels(regions.short$Region))
