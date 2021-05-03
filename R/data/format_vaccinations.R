## FORMATS VACCINATION DATA DIRECTLY FROM THE MODELLING CELL DRIVE
## This is a time consuming job (large data over a remote connection), so will only do this if appropriately named output files are not present
## or if a data overwrite flag has been set.

## ============ Where to find the data? ============
## -------------------------------------------------
## Either already specified in the variable vacc.loc
## or found as the latest entry into a default
## directory. If vacc.loc is NULL, look for a
## command-line argument.
require(tidyverse)
require(cubelyr)
require(lubridate)
if(!exists("vacc.loc")){ ## Set to default format for the filename
    input.loc <- "~/CoVID-19/Data streams/Vaccine line list"
    ## List the possible files in the directory
    vacc.loc <- file.info(file.path(input.loc,
                                    list.files(path = input.loc,
                                               pattern=glob2rx(paste("202", "immunisations SPIM", "zip", sep = "*")))
                                    )
                          )
    ## Has a particular date's data been specified
    rnames <- rownames(vacc.loc)
    if(exists("str.date.vacc")){
        input.loc <- rnames[grepl(str.date.vacc, rnames)]
    } else {
        ## Pick up the most recently added
        input.loc <- rnames[which.max(vacc.loc$ctime)]
        ## input.loc <- "~/CoVID-19/Data streams/Vaccine line list/20210212 immunisations SPIM.csv"
        str.date.vacc <- strapplyc(input.loc, "[0-9]{8,}", simplify = TRUE)
    }
} else {
    if(is.null(vacc.loc)){
        input.loc <- commandArgs(trailingOnly = TRUE)[1]
    } else {
        if(startsWith(vacc.loc, "/")) input.loc <- prev.loc
        else input.loc <- build.data.filepath(subdir = "raw", "vaccination", prev.loc)
    }
}

## Where will outputs be stored, to avoid repeat accessing of the remote COVID directory
vacc.rdata <- build.data.filepath(file.path("RTM_format", region.type, "vaccination"), region.type, "vacc", str.date.vacc, ".RData")

## Define an age-grouping
if(!exists("age.agg")){
    age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
    age.labs <- c("<1yr", "1-4", "5-14", "15-24", "25-44", "45-64", "65-74", "75+")
}
nA <- length(age.labs)

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

## Function to repeat the last row of a dataset until it pads out to the end of the series
pad.rows.at.end <- function(df, nrows.full){
    if(nrow(df) < nrows.full){
        nadds <- nrows.full - nrow(df)
        nold <- nrow(df)
        df[nold + (1:nadds), ] <- df[nold, ]
        df$sdate[nold + (1:nadds)] <- df$sdate[nold] + 1:nadds
    }
    df
}

## Function to format aggregated counts of vaccinations into the cross-tabulated datasets required for the rtm.
fn.region.crosstab <- function(dat, reg_r, dose_d, ndays = ndays){
    dat %>%
        filter(region == reg_r, dose == dose_d) %>%
        arrange(sdate) %>%
        group_by(age.grp, pop) %>%
        summarise(sdate, n, n.cum = cumsum(n)) %>%
        mutate(pop = max(pop, n.cum + 1)) %>%
        mutate(value = n / (pop - dplyr::lag(n.cum))) %>% ## Calculating the fraction of the denominator population still at risk.
        replace_na(list(value = 0)) %>%
        ungroup() %>%
        mutate(value = 2 * (1 - sqrt(1 - value))) %>% ## This line transforms the number of events until a final ok, and then 
        pivot_wider(id_cols = sdate,
                    names_from = age.grp,
                    values_from = value) %>%
        pad.rows.at.end(ndays)
    }

## Function to handle the unzipping of a large file - apparently R's unzip function is unreliable.
decompress_file <- function(directory, file, .file_cache = FALSE) {

    if (.file_cache == TRUE) {
       print("decompression skipped")
    } else {

        ## Set working directory for decompression
        ## simplifies unzip directory location behavior
        wd <- getwd()
        setwd(directory)

        ## Run decompression
        decompression <-
            system2("unzip",
                    args = c("-o", # include override flag
                             gsub(" ", "\\\\ ", file)),
                    stdout = TRUE)
        
        ## uncomment to delete archive once decompressed
        file.remove(file) 

        ## Reset working directory
        setwd(wd); rm(wd)

        ## Test for success criteria
        ## change the search depending on 
        ## your implementation
        if (grepl("Warning message", tail(decompression, 1))) {
            print(decompression)
        }

        return(gsub(".zip", ".csv", file.path(directory, file)))
        
    }
}

## Substitute this into the names of the intended data file names
vac1.files <- gsub("date.vacc", str.date.vacc, vac1.files, fixed = TRUE)
vacn.files <- gsub("date.vacc", str.date.vacc, vacn.files, fixed = TRUE)

## If these files exist and we don't want to overwrite them: do nothing
if(vac.overwrite || !all(file.exists(c(vac1.files, vacn.files)))){

    ## Extract file from archive
    if(!file.exists(gsub(".zip", ".csv", file.path("data", basename(input.loc))))){
        file.copy(input.loc, file.path("data", basename(input.loc)))
        input.loc <- decompress_file("data", basename(input.loc))
    } else input.loc <- gsub(".zip", ".csv", file.path("data", basename(input.loc)))
    
    if(str.date.vacc != strapplyc(input.loc, "[0-9]{8,}", simplify = TRUE)) stop("Specified date.vacc does not match the most recent prevalence data file.")
    
    ## Map our names for columns (LHS) to data column names (RHS)
    possible.col.names <- list(
        age = "age",
        region = c("region_of_residence", "Region_of_residence"),
        sdate = "vaccination_date",
        type = "product_display_type",
        dose = "string_dose_number",
        ltla_code = "ltla_code")
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
    vac.nhs.regions <- function(x) {
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
    
    if (region.type == "NHS") get.region <- vac.nhs.regions
    if (region.type == "ONS") get.region <- function(x) {
        y <- ons.region(x)
        x %>% mutate(region = y)
    }
    
    if(!exists("vac1.files")){
        vac1.files <- build.data.filepath("RTM_format/vaccination",
                                          str.date.vacc,
                                          "_1stvaccinations_",
                                          regions,
                                          ".txt")
        vacn.files <- build.data.filepath("RTM_format/vaccination",
                                          str.date.vacc,
                                          "_nthvaccinations_",
                                          regions,
                                          ".txt")
    }
    
    ## Which columns are we interested in?
    vacc.col.args <- list()
    vacc.col.args[[col.names[["type"]]]] <- col_character()
    vacc.col.args[[col.names[["age"]]]] <- col_double()
    vacc.col.args[[col.names[["region"]]]] <- col_character()
    vacc.col.args[[col.names[["dose"]]]] <- col_character()
    vacc.col.args[[col.names[["sdate"]]]] <- col_date(format="")
    vacc.col.args[[col.names[["ltla_code"]]]] <- col_character()
    vacc.cols <- do.call(cols, vacc.col.args)

    ## Reading in the data ##
    print(paste("Reading from", input.loc))
    ## strPos <- c("+", "Positive", "positive")
    vacc.dat <- read_csv(input.loc,
                         col_types = vacc.cols) %>%
        select(!!(names(vacc.cols$cols)))
    cat("Got here2\n")
    vacc.dat <- vacc.dat %>%
        rename(!!!col.names)
    cat("Got here 2a\n")
    vacc.dat <- vacc.dat %>%
	mutate(sdate = fuzzy_date_parse(sdate),
               age.grp = cut(age, age.agg, age.labs, right = FALSE, ordered_result = T)
               )
    cat("Got here3\n")
    ## Remove rows with missing key information
    n.vaccs <- nrow(vacc.dat)
    vacc.dat <- vacc.dat %>% filter(region != "Unknown",
                                    !is.na(type),
                                    !is.na(sdate),
                                    !is.na(dose),
                                    !is.na(age.grp))
    n.vaccs.complete <- nrow(vacc.dat)
    cat("Got here4\n")
    ## r.even <- function(vaccs, len) rmultinom(1, vaccs, rep(1, len))
    
    ## ## ====== DATA ON FIRST VACCINATION
    
    ## vaccs <- rbind(r.even(269, 20), r.even(11760, 20), r.even(9327, 10), r.even(5746, 5), r.even(5681, 5), r.even(4153, 5), r.even(1361, 5), r.even(583, 5), r.even(386, 5), r.even(14424, 1))
    ## vaccs <- cbind(vaccs, rbind(r.even(1352, 20), r.even(45379, 20), r.even(36858, 10), r.even(23132, 5), r.even(23686, 5), r.even(17056, 5), r.even(6733, 5), r.even(4228, 5), r.even(7299, 5), r.even(395103, 1)))
    ## vaccs <- cbind(vaccs, rbind(r.even(2202, 20), r.even(75231, 20), r.even(60873, 10), r.even(37872, 5), r.even(38787, 5), r.even(27809, 5), r.even(11101, 5), r.even(7071, 5), r.even(13338, 5), r.even(543397, 1)))
    ## vaccs <- cbind(vaccs, rbind(r.even(3563, 20), r.even(124955, 20), r.even(97479, 10), r.even(59808, 5), r.even(60157, 5), r.even(43763, 5), r.even(18352, 5), r.even(13783, 5), r.even(32634, 5), r.even(670040, 1)))
    ## dimnames(vaccs) <- list(ages = 0:80, date = lubridate::as_date(c("20201213","20201220","20201227","20210103")))
    
    ## vaccs <- vaccs %>% as.tbl_cube(met_name = "jabs") %>% as_tibble() %>% mutate(date = lubridate::as_date(date), Age.Grp = cut(ages, breaks = c(0, 1, 5, 15, 25, 45, 65, 75, Inf), right = FALSE, ordered_result = T)) %>% group_by(Age.Grp, date) %>% summarise(jabs = sum(jabs))
    
    ## probs <- c(8942, 4969, 13188, 5172, 7416, 7709, 6218)
    ## probs <- rbind(probs, c(70593, 59141, 112177, 78726, 71914, 92577, 75194), c(86555, 80364, 152828, 145645, 106357, 143259, 101894), c(112616, 115341, 211936, 201823, 144815, 195125, 141767))
    ## dimnames(probs) <- list(date = lubridate::as_date(c("20201213","20201220","20201227","20210103")), region = regions)
    ## probs <- apply(probs, 1, function(x) x / sum(x))
    ## probs <- probs %>% as.tbl_cube(met_name = "propn") %>% as_tibble() %>% mutate(date = lubridate::as_date(date))
    
    ## merged <- vaccs %>% left_join(probs)
    ## merged_wide <- pivot_wider(merged, id_cols = 1:3, names_from = region, values_from = propn)
    
    ## tmp.samples <- NULL
    ## for(i in 1:nrow(merged_wide))
    ##     tmp.samples <- cbind(tmp.samples, rmultinom(1, merged_wide[i, ]$jabs, c(merged_wide[i, ]$East_of_England, merged_wide[i, ]$London, merged_wide[i, ]$Midlands, merged_wide[i, ]$North_East, merged_wide[i, ]$North_West, merged_wide[i, ]$South_East, merged_wide[i, ]$South_West)))
    ## rownames(tmp.samples) <- colnames(merged_wide)[-(1:3)]
    ## sampled.jabs <- merged_wide[, 1:2]
    ## tmp.samples <- as.data.frame(t(tmp.samples))
    ## sampled.jabs <- bind_cols(sampled.jabs, tmp.samples)
    
    ## sampled.jabs <- sampled.jabs %>% pivot_longer(cols = -(1:2), names_to = "Region", values_to = "Jabs") %>% right_join(expand.grid(date = lubridate::ymd("20200217"):lubridate::ymd("20210117"), Region = colnames(merged_wide)[-(1:3)], Age.Grp = unique(merged_wide$Age.Grp)) %>% as.data.frame %>% mutate(`date` = lubridate::as_date(`date`))) %>% replace_na(list(Jabs = 0)) %>% arrange(date)
    
    names(vac1.files) <- names(vacn.files) <- regions
    vac.dates <- as_date(earliest.date:(max(vacc.dat$sdate) + vacc.lag))
    jab.dat <- vacc.dat %>%
        group_by(sdate, ltla_code, region, age.grp, dose) %>%
        summarise(n = n(), pPfizer = sum(type == "Pfizer") / n()) %>%
        ungroup() %>%
        mutate(sdate = sdate + vacc.lag) %>%
        get.region() %>%
        select(-ltla_code)
    rm(vacc.dat)
    ## if(region.type == "ONS")
    jab.dat <- jab.dat %>%
        group_by(sdate, region, age.grp, dose) %>%
        summarise(pPfizer = weighted.mean(pPfizer, n), n = sum(n))
    jab.dat <- jab.dat %>%
        right_join(expand_grid(sdate = vac.dates,
                               region = regions,
                               age.grp = factor(age.labs, levels = age.labs, ordered = TRUE),
                               dose = c("First", "Second")),
                   by = c("sdate", "region", "age.grp", "dose")
                   ) %>%
        replace_na(list(n = 0, pPfizer = 1)) %>% ## What to fill in the columns where there are no data
        arrange(sdate) %>%
        left_join(pdf.all) %>%
        rename(pop = count)
    cat("Got here5\n")
    ## Will need to extract some design matrices for vaccine efficacy also
    v1.design <- NULL
    vn.design <- NULL

    ## Following code will extrapolate the vaccination programme out into the future.
    source(file.path(proj.dir, "R", "data", "augment_vaccinations.R"))

    for(reg in regions){
        region.dat <- jab.dat %>% ungroup() %>% fn.region.crosstab(reg, "First", ndays = max(jab.dat$sdate) - ymd("20200216"))
        cat("First stopifnot\n")
        stopifnot(!is.na(region.dat))
        stopifnot(all(region.dat[, -1] >= 0))
        stopifnot(all(region.dat[, -1] <= 2))
        tmpFile <- vac1.files[reg]
        dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)
        region.dat %>%
            write_tsv(tmpFile,
                      col_names = FALSE)

        if(vac.design == "cumulative"){
            tmp.design <- as.vector(
                as.matrix(
                    pivot_wider(jab.dat %>%
                                filter(region == reg, dose == "First", sdate %in% vac.dates) %>%
                                arrange(sdate) %>%
                                group_by(age.grp, pop) %>%
                                summarise(pPfizer.wt = ifelse(cumsum(n) == 0, 1, cumsum(n * pPfizer) / cumsum(n)), sdate = sdate),
                                id_cols = age.grp, names_from = sdate, values_from = pPfizer.wt) %>%
                    ungroup() %>%
                    select(-age.grp)
                )
            )
        } else
            tmp.design <- as.vector(as.matrix(pivot_wider(jab.dat %>% filter(region == reg, dose == "First", sdate %in% vac.dates), id_cols = age.grp, names_from = sdate, values_from = pPfizer) %>% select(-age.grp)))
        
        v1.design <- rbind(v1.design, cbind(tmp.design, 1 - tmp.design))

        region.dat <- jab.dat %>%
            filter(region == reg, dose == "Second") %>%
            left_join(jab.dat %>%
                      filter(region == reg, dose == "First") %>%
                      arrange(sdate) %>%
                      group_by(region, age.grp) %>%
                      summarise(sdate = sdate, sum.first = cumsum(n))
                      ) %>%
            arrange(sdate) %>%
            group_by(age.grp) %>%
            summarise(sdate = sdate,
                      value = n / (dplyr::lag(sum.first) - dplyr::lag(cumsum(n)))
                      ) %>%
            replace_na(list(value = 0)) %>%
            ungroup() %>%
            mutate(value = 2 * (1 - sqrt(1 - value))) %>%   ## This line transforms the data from a fraction to a rate
            pivot_wider(id_cols = sdate,
                        names_from = age.grp,
                        values_from = value) %>%
            pad.rows.at.end(ndays)
        cat("Second stopifnot\n")
        stopifnot(!is.na(region.dat))
        stopifnot(all(region.dat[, -1] >= 0))
        stopifnot(all(region.dat[, -1] <= 2))
        tmpFile <- vacn.files[reg]
        dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)
        region.dat %>%
            write_tsv(tmpFile,
                      col_names = FALSE)

        if(vac.design == "cumulative"){
            tmp.design <- as.vector(
                as.matrix(
                    pivot_wider(jab.dat %>%
                                filter(region == reg, dose == "Second", sdate %in% vac.dates) %>%
                                arrange(sdate) %>%
                                group_by(age.grp, pop) %>%
                                summarise(pPfizer.wt = ifelse(cumsum(n) == 0, 1, cumsum(n * pPfizer) / cumsum(n)), sdate = sdate),
                                id_cols = age.grp, names_from = sdate, values_from = pPfizer.wt) %>%
                    ungroup() %>%
                    select(-age.grp)
                )
            )
        } else 
            tmp.design <- as.vector(as.matrix(pivot_wider(jab.dat %>% filter(region == reg, dose == "Second", sdate %in% vac.dates), id_cols = age.grp, names_from = sdate, values_from = pPfizer) %>% select(-age.grp)))
        
        vn.design <- rbind(vn.design, cbind(tmp.design, 1 - tmp.design))
        
    }
    
    ## ## ====== DATA ON SECOND VACCINATION
    
    ## vaccs <- rbind(r.even(0, 20), r.even(0, 20), r.even(0, 10), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 1))
    ## vaccs <- cbind(vaccs, rbind(r.even(0, 20), r.even(0, 20), r.even(0, 10), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 1)))
    ## vaccs <- cbind(vaccs, rbind(r.even(0, 20), r.even(0, 20), r.even(0, 10), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 1)))
    ## vaccs <- cbind(vaccs, rbind(r.even(81, 20), r.even(3700, 20), r.even(3125, 10), r.even(1951, 5), r.even(1876, 5), r.even(1311, 5), r.even(476, 5), r.even(240, 5), r.even(153, 5), r.even(6156, 1)))
    ## dimnames(vaccs) <- list(ages = 0:80, date = lubridate::as_date(c("20201213","20201220","20201227","20210103")))
    
    ## vaccs <- vaccs %>% as.tbl_cube(met_name = "jabs") %>% as_tibble() %>% mutate(date = lubridate::as_date(date), Age.Grp = cut(ages, breaks = c(0, 1, 5, 15, 25, 45, 65, 75, Inf), right = FALSE, ordered_result = T)) %>% group_by(Age.Grp, date) %>% summarise(jabs = sum(jabs))
    
    ## probs <- rep(1, 7)
    ## probs <- rbind(probs, rep(1, 7), rep(1,7), c(3710, 1102, 4165, 1102, 4501, 2408, 2041))
    ## dimnames(probs) <- list(date = lubridate::as_date(c("20201213","20201220","20201227","20210103")), region = regions)
    ## probs <- apply(probs, 1, function(x) x / sum(x))
    ## probs <- probs %>% as.tbl_cube(met_name = "propn") %>% as_tibble() %>% mutate(date = lubridate::as_date(date))
    
    ## merged <- vaccs %>% left_join(probs)
    ## merged_wide <- pivot_wider(merged, id_cols = 1:3, names_from = region, values_from = propn)
    
    ## tmp.samples <- NULL
    ## for(i in 1:nrow(merged_wide))
    ##     tmp.samples <- cbind(tmp.samples, rmultinom(1, merged_wide[i, ]$jabs, c(merged_wide[i, ]$East_of_England, merged_wide[i, ]$London, merged_wide[i, ]$Midlands, merged_wide[i, ]$North_East, merged_wide[i, ]$North_West, merged_wide[i, ]$South_East, merged_wide[i, ]$South_West)))
    ## rownames(tmp.samples) <- colnames(merged_wide)[-(1:3)]
    ## sampled.jabs <- merged_wide[, 1:2]
    ## tmp.samples <- as.data.frame(t(tmp.samples))
    ## sampled.jabs <- bind_cols(sampled.jabs, tmp.samples)
    
    ## sampled.jabs <- sampled.jabs %>% pivot_longer(cols = -(1:2), names_to = "Region", values_to = "Jabs") %>% right_join(expand.grid(date = lubridate::ymd("20200217"):lubridate::ymd("20210117"), Region = colnames(merged_wide)[-(1:3)], Age.Grp = unique(merged_wide$Age.Grp)) %>% as.data.frame %>% mutate(`date` = lubridate::as_date(`date`))) %>% replace_na(list(Jabs = 0)) %>% arrange(date)
    
    ## vacc2.files <- paste0(file.path(getwd(), "data", "RTM_format", "NHS", "vaccination", "dummy_dose2_"), regions, ".txt")
    ## names(vacc2.files) <- regions
    
    ## for(reg in regions){
    ##     region.dat <- pivot_wider(sampled.jabs %>% filter(Region == reg),
    ##                               id_cols = 2,
    ##                               names_from = Age.Grp,
    ##                               values_from = Jabs)
    ##     tmpFile <- vacc2.files[reg]
    ##     dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)
        
    ##     region.dat %>%
    ##         write_tsv(tmpFile,
    ##                   col_names = FALSE)
    ## }

    save(vac.dates, v1.design, vn.design, jab.dat, file = vacc.rdata)
    ## file.remove(file.path("data", basename(input.loc)))
    
} else {

    load(vacc.rdata)
    
}
