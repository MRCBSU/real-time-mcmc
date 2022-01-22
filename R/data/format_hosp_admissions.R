suppressMessages(require(lubridate))
suppressMessages(require(tidyverse))
suppressMessages(require(readxl))
suppressMessages(require(gsubfn))

#########################################################
## Inputs that should (or may) change on a weekly basis
#########################################################
## Where to find the data, if NULL use command line argument
if(!exists("adm_seb.loc") | !exists("adm_sus.loc")){ ## Set to default format for the filename
    #input.loc <- build.data.filepath(subdir = "raw", "serology")
    input.loc <- "data/raw/admissions"
    ## List the possible files in the directory
    adm_seb.loc <- file.info(file.path(input.loc, list.files(path=input.loc, pattern=glob2rx("sitrep*rds"))))
    if (region.type == "NHS") {
        adm_sus.loc <- file.info(file.path(input.loc, list.files(path=input.loc, pattern=glob2rx("nhs*csv"))))
    } else {
        adm_sus.loc <- file.info(file.path(input.loc, list.files(path=input.loc, pattern=glob2rx("ltla*csv"))))
    }
    ## Pick the most recently added
    input_seb.loc <- rownames(adm_seb.loc)[which.max(adm_seb.loc$mtime)]
    input_sus.loc <- rownames(adm_sus.loc)[which.max(adm_sus.loc$mtime)]
} else {
    if(is.null(adm_sus.loc) | is.null(adm_seb.loc)){
        input_seb.loc <- commandArgs(trailingOnly = TRUE)[1]
        input_sus_ons.loc <- commandArgs(trailingOnly = TRUE)[2]
    } else {
        if(startsWith(adm_seb.loc, "/")) input_seb.loc <- adm_seb.loc
        else input_seb.loc <- build.data.filepath(subdir = "raw", "admissions", adm_seb.loc)
        if(startsWith(adm_sus.loc, "/")) input_sus.loc <- adm_sus.loc
        else input_sus.loc <- build.data.filepath(subdir = "raw", "admissions", adm_sus.loc)
    }
}

## Are we using SUS data that has been pre-processed once before?
if(sus_seb_combination == 3) {
    old_input_adm <- "data/previous_run_input/"
    if(region.type == "NHS") {
        old_adm.loc <- file.path(old_input_adm, "NHS", "admissions", preprocessed_sus_names)
        old_adm_csv.loc <- file.path(old_input_adm, "NHS", "admissions", preprocessed_sus_csv_name)
    } else {
        old_adm.loc <- file.path(old_input_adm, "ONS", "admissions", preprocessed_sus_names)
        old_adm_csv.loc <- file.path(old_input_adm, "ONS", "admissions", preprocessed_sus_csv_name)
    }
    names(old_adm.loc) <- names(preprocessed_sus_names)
}

## What is the date of publication of these data? If not specified, try to extract from filename
#! Note get peter to output sus data in this format going forwards
if(!exists("date.adm_sus")){
    fl.name <- basename(input.loc)
    date.adm_sus.str <- str_match(rownames(adm_sus.loc)[which.max(adm_sus.loc$mtime)],"([0-9]+)\\.csv$")[,2]
    date.adm_sus <- ymd(date.adm_sus.str)
}
if(!exists("date.adm_seb")){
    fl.name <- basename(input.loc)
    date.adm_seb.str <- str_match(rownames(adm_seb.loc)[which.max(adm_seb.loc$mtime)],"([0-9]+)\\.rds$")[,2]
    date.adm_seb <- ymd(date.adm_seb.str)
}

# Set the dates to start and end the different sections of data (sus vs sus + seb vs seb) )
earliest.date <- start.date
if(sus_seb_combination %in% c(1, 3)) {
    latest_sus.date <- adm_sus.end.date
    earliest_seb.date <- adm_sus.end.date + 1
}
if(exists("adm.end.date")){
    latest.date <- adm.end.date
} else {
    latest.date <- date.adm_seb - adm_seb.strip_days
}

## Define an age-grouping (the final age groupings note these can be modified in config.r this is just a default)
if(!exists("age_adm.agg")){
    age_adm.agg <- c(0, 25, 45, 65, 75, Inf)
    ## age_adm.oldlabs <- c("0_25", "25_45", "45_65", "65_75", "75", "")
    age_adm.labs <- c("0-25", "25-45","45-65", "65-75", "75+", "NA") ## "All ages"
} ## else age_adm.oldlabs <- "Needs defining"
nA_adm <- length(age_adm.labs)

## Read in sebs data if it is needed
## Note data is stored in an rds so we must initially read in all data
if(sus_seb_combination %in% c(1, 2, 3)) {

    # Only list the values with age data for admissions and diagnoses
    possible.col.names.seb <- list(
        region = c("region", "Region", "Geography", "NHS_Region"),
        date = c("DateVal", "date"),
        trust_code = c("org_code"),
        trust_name = c("org_name"),
        # admissions = c("n_patients_admitted"),
        admissions_ages_0_5 = c("n_patients_admitted_age_0_5"),
        admissions_ages_6_17 = c("n_patients_admitted_age_6_17"),
        admissions_ages_18_24 = c("n_patients_admitted_age_18_24"),
        admissions_ages_25_34 = c("n_patients_admitted_age_25_34"),
        admissions_ages_35_44 = c("n_patients_admitted_age_35_44"),
        admissions_ages_45_54 = c("n_patients_admitted_age_45_54"),
        admissions_ages_55_64 = c("n_patients_admitted_age_55_64"),
        admissions_ages_65_74 = c("n_patients_admitted_age_65_74"),
        admissions_ages_75_84 = c("n_patients_admitted_age_75_84"),
        admissions_ages_85 = c("n_patients_admitted_age_85"),
        # diagnoses = c("n_inpatients_diagnosed"),
        diagnoses_ages_0_5 = c("n_inpatients_diagnosed_age_0_5"),
        diagnoses_ages_6_17 = c("n_inpatients_diagnosed_age_6_17"),
        diagnoses_ages_18_24 = c("n_inpatients_diagnosed_age_18_24"),
        diagnoses_ages_25_34 = c("n_inpatients_diagnosed_age_25_34"),
        diagnoses_ages_35_44 = c("n_inpatients_diagnosed_age_35_44"),
        diagnoses_ages_45_54 = c("n_inpatients_diagnosed_age_45_54"),
        diagnoses_ages_55_64 = c("n_inpatients_diagnosed_age_55_64"),
        diagnoses_ages_65_74 = c("n_inpatients_diagnosed_age_65_74"),
        diagnoses_ages_75_84 = c("n_inpatients_diagnosed_age_75_84"),
        diagnoses_ages_85 = c("n_inpatients_diagnosed_age_85")
    )

    # Ensure all the useful pieces of data have been grabbed
    adm.dat.seb <- read_rds(input_seb.loc)
    input.col.names.seb <- suppressMessages(names(adm.dat.seb))
    is.valid.col.name.seb <- function(name) {name %in% input.col.names.seb}
    first.valid.col.name.seb <- function(names) {first.where.true(names, is.valid.col.name.seb)}
    col.names.seb <- lapply(possible.col.names.seb, first.valid.col.name.seb)
    invalid.col.names <- sapply(col.names.seb, is.null)
    if (any(invalid.col.names)) {
        names.invalid.cols <- paste0(names(possible.col.names.seb)[invalid.col.names], collapse = ", ")
        stop(paste("No valid column name for:", names.invalid.cols))
    }
}

# Select only the columns we want to read in from the sus data (if we need the sus data)
if(sus_seb_combination %in% c(0,1)) {
    possible.col.names.sus <- list(
        date = c("DateVal", "date"),
        ages = c("ageGrpRTM", "ages"),
        nosocomial = c("group")
    )


    input.col.names.sus <- suppressMessages(names(read_csv(input_sus.loc, n_max = 0, na = "")))
    is.valid.col.name.sus <- function(name) {name %in% input.col.names.sus}
    first.valid.col.name.sus <- function(names) {first.where.true(names, is.valid.col.name.sus)}
    col.names.sus <- lapply(possible.col.names.sus, first.valid.col.name.sus)
    invalid.col.names <- sapply(col.names.sus, is.null)
    if (any(invalid.col.names)) {
        names.invalid.cols <- paste0(names(possible.col.names.sus)[invalid.col.names], collapse = ", ")
        stop(paste("No valid column name for:", names.invalid.cols))
    }
}

# define get region functions (note modified to take in the column rather than the df)
get.region <- function(x) str_replace_all(x, " ", "_")

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

if(!exists("admsam.files")){
    if(!exists("adm.end.date")){
        date.adm.str <- lubridate::as_date(ifelse(sus_seb_combination > 0,
                                                  date.adm_seb - adm_seb.strip_days,
                                                  date.adm_sus - adm_sus.strip_days))
    } else date.adm.str <- adm.end.date
    admsam.files <- paste0(data.dirs["adm"], "/", date.adm.str, "_", regions, "_", nA_adm, "ag_counts.txt")
}

## Construct the sus data into a useful format if necessary (date, age, region, admissions)
if(sus_seb_combination %in% c(0, 1)){

    ## Define the correct column types
    adm.date.fmt <- "%y-%b-%d"
    adm.col.args <- list()
    adm.col.args[[col.names.sus[["date"]]]] <- col_date(adm.date.fmt)
    adm.col.args[[col.names.sus[["ages"]]]] <- col_character()
    adm.col.args[[col.names.sus[["nosocomial"]]]] <- col_character()
    adm.cols <- do.call(cols_only, adm.col.args)

    # define a function to get the col_types
    #! Note issue with default being integer (currently necessary as ons +nhs data uses many regions, ons is worst as would have to define all ltla regions as columns manually)
    get_col_type_sus <- function(x) {
        ifelse(
            identical(x, col_character()), "c",
            ifelse(identical(x, col_double()), "d",
            ifelse(identical(x, col_integer()), "i",
            ifelse(identical(x, col_date(adm.date.fmt)), "D",
            "i"
        ))))
    }
    # get the coltypes and name all columns with df colnames
    fields <- sapply(adm.col.args[input.col.names.sus], get_col_type_sus)
    names(fields) <- ifelse(is.na(names(fields)), input.col.names.sus, names(fields))

    # Read the data from the sus file and rearrange to summarise across ages
    print(paste0("Reading in data from ", input.loc))

    adm.dat.sus <- read_csv(input_sus.loc, col_types = fields, na = "")  %>%
          rename(!!!col.names.sus)  %>%
        pivot_longer(where(is.numeric), names_to = "region", values_to = "admissions") %>%
        ungroup()  %>%
        group_by(date, ages, region)  %>%
        summarise(across(where(is.numeric), ~sum(., na.rm = T)))  %>%
        pivot_wider(id_cols = c(date, region), names_from = ages, values_from = admissions)

    # Combine ages as defined in the config file
    for(i in seq_along(summarise_classes_sus)) {
        adm.dat.sus <- adm.dat.sus %>%
                group_by(date, region)  %>%
                mutate(!!names(summarise_classes_sus)[[i]] := sum(!!!syms(summarise_classes_sus[[i]]), na.rm = T), .keep = "unused")
    }

    # Give ages correct names (match with both dfs)
    adm.dat.sus <- adm.dat.sus  %>%
        pivot_longer(where(is.numeric), names_to = "ages", values_to = "admissions")  %>%
        mutate(ages = factor(ages, levels = age_adm_sus.oldlabs))

    levels(adm.dat.sus$ages) <- age_adm.labs

    # Relabel problem regions as defined in config
    #! Note issue that it assumes consistency of region names within a df
    if(region.type == "NHS") {

        adm.dat.sus <- adm.dat.sus %>%
            mutate(region = as.factor(get.region(region)))

        if(length(relevel_sus_locs_nhs) != 0) {
            adm.dat.sus  <- adm.dat.sus  %>%
                mutate(region = fct_recode(region, !!!relevel_sus_locs_nhs))
        }
    } else {

        # Read in linker with linking column and new region column defined in config.r
        possible.col.names.sus_link <- list(
            link = c(adm.sus.geog_link),
            region = c(adm.sus.region_col)
        )

        # Check the columns are correct
        input.col.names.sus_link <- suppressMessages(names(read_csv(paste(input.loc,adm.sus.geog_link.loc, sep = "/"), n_max = 0)))
        is.valid.col.name.sus_link <- function(name) {name %in% input.col.names.sus_link}
        first.valid.col.name.sus_link <- function(names) {first.where.true(names, is.valid.col.name.sus_link)}
        col.names.sus_link <- lapply(possible.col.names.sus_link, first.valid.col.name.sus_link)
        invalid.col.names <- sapply(col.names.sus_link, is.null)
        if (any(invalid.col.names)) {
            names.invalid.cols <- paste0(names(possible.col.names.sus_link)[invalid.col.names], collapse = ", ")
            stop(paste("No valid column name for:", names.invalid.cols))
        }

        # define column types
        adm.date.fmt <- "%y-%b-%d"
        adm.col.args <- list()
        adm.col.args[[col.names.sus_link[["link"]]]] <- col_character()
        adm.col.args[[col.names.sus_link[["region"]]]] <- col_character()
        adm.cols <- do.call(cols_only, adm.col.args)

        # define column types and defaults as function
        get_col_type_sus_link <- function(x) {
            ifelse(
                identical(x, col_character()), "c",
                ifelse(identical(x, col_double()), "d",
                ifelse(identical(x, col_integer()), "i",
                ifelse(identical(x, col_date(adm.date.fmt)), "D",
                "_"
            ))))
        }

        # set coltypes and give coltypes correct colnames
        fields <- sapply(adm.col.args[input.col.names.sus_link], get_col_type_sus_link)
        names(fields) <- ifelse(is.na(names(fields)), input.col.names.sus_link, names(fields))

        # Read in geography linker
        sus_geog_link <- read_csv(paste(input.loc,adm.sus.geog_link.loc, sep = "/"), col_types = fields)  %>%
                            rename(!!!col.names.sus_link)

        # Join on linker and select new region column as only geography column
        adm.dat.sus <- sus_geog_link %>%
                        right_join(adm.dat.sus, by = c("link" = "region")) %>%
                        mutate(region = get.region(region))  %>%
                        select(-link)

        # Give ages correct names (match with both dfs)
        if(length(relevel_sus_locs_ons) != 0) {
            adm.dat.sus  <- adm.dat.sus  %>%
                mutate(region = fct_recode(region, !!!relevel_sus_locs_ons))
        }
    }

    # Sum all admissions by each group
    adm.dat.sus <- adm.dat.sus %>%
        group_by(region, date, ages)  %>%
        summarise(across(admissions, ~sum(., na.rm = T)))

}


## Construct sebs data into a useful format if necessary (date, age, region, admissions)
if(sus_seb_combination %in% c(1, 2, 3)) {
    # Rearrange sebs data to summarise across the ages (grab ages from the columns)
    adm.dat.seb <- adm.dat.seb  %>%
        select_if(!(names(.) %in% c("n_beds_mechanical_non_cov_ns", "n_beds_mechanical_suspected", "n_beds_noninvasive_non_cov_ns", "n_beds_noninvasive_suspected", "n_beds_oxygenation_non_cov_ns", "n_beds_oxygenation_suspected", "n_beds_other_3_non_cov_ns", "n_beds_other_3_suspected", "suspected_prev", "non_cov_nsid_prev", "non_cov_ns_prev", "n_staff_absent_caring", "n_staff_absent_positive", "n_care_home_admitted"))) %>% 
        unique()  %>%
        rename(!!!col.names.seb)  %>%
        select(!!!(unlist(names(col.names.seb)))) %>%
        pivot_longer(where(is.numeric), names_to = "ages", values_to = "admissions") %>%
        mutate(ages = str_match(ages, "_ages\\_(.+)$")[,2])  %>%
        group_by(ages, date, region, trust_code, trust_name)  %>%
        summarise(across(admissions, ~sum(., na.rm = T)))  %>%
        pivot_wider(id_cols = c(date, region, trust_code, trust_name), names_from = ages, values_from = admissions)

    # Combine ages as defined in the config file
    for(i in seq_along(summarise_classes_seb)) {
        adm.dat.seb <- adm.dat.seb %>%
                        group_by(date,region, trust_code, trust_name)  %>%
                        mutate(!!names(summarise_classes_seb)[[i]] := sum(!!!syms(summarise_classes_seb[[i]]), na.rm = T), .keep = "unused")
    }

    # Give ages the correct names (match with both DFs)
    adm.dat.seb <- adm.dat.seb  %>%
    pivot_longer(where(is.numeric), names_to = "ages", values_to = "admissions")  %>%
        mutate(ages = factor(ages, levels = age_adm_seb.oldlabs))

    levels(adm.dat.seb$ages) <- age_adm.labs

    # Relabel problem regions as defined in config
    #! Note issue that it assumes consistency of region names within a df
    if(region.type == "NHS") {

        adm.dat.seb <- adm.dat.seb %>%
            mutate(region = as.factor(get.region(region)))

        if(length(relevel_seb_locs_nhs) != 0) {
            adm.dat.seb  <- adm.dat.seb  %>%
                mutate(region = fct_recode(region, !!!relevel_seb_locs_nhs))
        }

    } else {

        # Read in linker with linking column and new region column defined in config.r
        possible.col.names.seb_link <- list(
            link = c(adm.seb.geog_link),
            region = c(adm.seb.region_col)
        )

        # Check the columns are correct
        input.col.names.seb_link <- suppressMessages(names(read_xlsx(paste(input.loc,adm.seb.geog_link.loc, sep = "/"), n_max = 0)))
        is.valid.col.name.seb_link <- function(name) {name %in% input.col.names.seb_link}
        first.valid.col.name.seb_link <- function(names) {first.where.true(names, is.valid.col.name.seb_link)}
        col.names.seb_link <- lapply(possible.col.names.seb_link, first.valid.col.name.seb_link)
        invalid.col.names <- sapply(col.names.seb_link, is.null)
        if (any(invalid.col.names)) {
            names.invalid.cols <- paste0(names(possible.col.names.seb_link)[invalid.col.names], collapse = ", ")
            stop(paste("No valid column name for:", names.invalid.cols))
        }

        # define column types
        adm.date.fmt <- "%y-%b-%d"
        adm.col.args <- list()
        adm.col.args[[col.names.seb_link[["link"]]]] <- col_character()
        adm.col.args[[col.names.seb_link[["region"]]]] <- col_character()
        adm.cols <- do.call(cols_only, adm.col.args)

        # define column types and defaults as function
        get_col_type_seb_link <- function(x) {
            ifelse(
                identical(x, col_character()), "text",
                ifelse(identical(x, col_double()), "numeric",
                ifelse(identical(x, col_integer()), "numeric",
                ifelse(identical(x, col_date(adm.date.fmt)), "date",
                "skip"
            ))))
        }

        # set coltypes and give coltypes correct colnames
        fields <- sapply(adm.col.args[input.col.names.seb_link], get_col_type_seb_link)
        names(fields) <- ifelse(is.na(names(fields)), input.col.names.seb_link, names(fields))

        # Read in geography linker
        seb_geog_link <- read_xlsx(paste(input.loc,adm.seb.geog_link.loc, sep = "/"), col_types = fields)  %>%
                            rename(!!!col.names.seb_link)

        # Join on linker and select new region column as only geography column
        adm.dat.seb <- seb_geog_link %>%
                        right_join(adm.dat.seb %>% ungroup() %>% select(-region), by = c("link" = "trust_code")) %>%
                        filter(!is.na(region))  %>%
                        mutate(region = get.region(region))  %>%
                        select(-link)

        # Give ages correct names (match with both dfs)
        if(length(relevel_seb_locs_ons) != 0) {
            adm.dat.seb  <- adm.dat.seb  %>%
                mutate(region = fct_recode(region, !!!relevel_seb_locs_ons))
        }
    }

    # Sum all admissions by each group
    adm.dat.seb <- adm.dat.seb  %>%
        group_by(region, date, ages)  %>%
        summarise(across(admissions, ~sum(., na.rm = T)))  %>%
        mutate(date = date - seb_report_delay)

}

# Combine sebs data and sus data into the same df and select only correct dates
if(sus_seb_combination == 0) {
    adm.sam <- adm.dat.sus %>%
                filter(date >= earliest.date)  %>%
                filter(date <= date.adm_sus - adm_sus.strip_days)  %>%
                filter(date <= latest.date)
} else if(sus_seb_combination == 1) {
    adm.sam <- adm.dat.sus %>%
                filter(date >= earliest.date)  %>%
                filter(date <= date.adm_sus - adm_sus.strip_days)  %>%
                filter(date <= latest_sus.date)  %>%
                bind_rows(adm.dat.seb  %>%
                            filter(date >= earliest_seb.date)  %>%
                            filter(date <= date.adm_seb - adm_seb.strip_days)  %>%
                            filter(date <= latest.date)
                )
} else if(sus_seb_combination == 2){
    adm.sam <- adm.dat.seb  %>%
        filter(date >= earliest.date)  %>%
        filter(date <= date.adm_seb - adm_seb.strip_days)  %>%
        filter(date <= latest.date)
} else if(sus_seb_combination == 3) {
    adm.sam <- adm.dat.seb %>%
        filter(date >= earliest_seb.dat) %>%
        filter(date <= date.adm_seb - adm_seb.strip_days) %>%
        filter(date <= latest.date)
}

## Write rtm data outputs to file
names(admsam.files) <- regions

## Create missing directories
walk(dirname(admsam.files), ~dir.create(., showWarnings = F, recursive = T))

## Pad data with some zeros to take it back to day 1.
adm.sam <- expand_grid(region = regions, date = lubridate::as_date(start.date:max(adm.sam$date)), ages = unique(adm.sam$ages)) %>%
    mutate(admissions = 0) %>%
    bind_rows(adm.sam) %>%
    group_by(region, date, ages) %>%
    summarise(admissions = max(admissions)) %>%
    ungroup()

# Loop over regions and save the data for them
for(reg in regions) {
    region.sam <- pivot_wider(adm.sam %>%
                              filter(region == reg),
                              id_cols = c(date, region),
                              names_from = ages,
                              values_from = admissions
                              ) %>%
        ungroup() %>%
        select(-region)

    if(sus_seb_combination == 3) {
        tmp_sus <- read_tsv(old_adm.loc[reg], col_names = F)
        colnames(tmp_sus) <- colnames(region.sam)
        
        tmp_sus <- tmp_sus %>%
            mutate(date = as_date(date)) %>%
            filter(date <= ymd(latest_sus.date))
        
        region.sam <- tmp_sus %>%
            bind_rows(
                region.sam %>% 
                filter(date >= latest_sus.date)
            )
        
    }

    tmpFile <- admsam.files[reg]

    print(paste("Writing to", tmpFile))

    region.sam %>%
        write_tsv(
            tmpFile,
            col_names = F
        )
}

# Create missing directory
if(!file.exists(out.dir)) dir.create(out.dir, recursive = T)

if(sus_seb_combination == 3) {
    tmp_sus <- read_csv(old_adm_csv.loc) %>%
        mutate(date = as_date(date)) %>%
        filter(date <= ymd(latest_sus.date))

    adm.sam <- tmp_sus %>%
            bind_rows(
                adm.sam %>% 
                    filter(date >= latest_sus.date)
            ) %>%
            arrange(region, date, ages, admissions)
}

# Save the data
write_csv(adm.sam, file.path(out.dir, "admissions_data.csv"))


## Create a quick plot of the data
require(ggplot2)
require(viridis)
require(hrbrthemes)
require(plotly)
require(htmlwidgets)

adm.plot <- adm.sam %>%
    filter(!is.na(region)) %>%
    group_by(date, region) %>%
    summarise(across(admissions, .fn = sum)) %>%
    mutate(text = paste("Region:", region, "\nAdmssions:", admissions, "\nDate:", date))

ggp <- ggplot(adm.plot, aes(x = date, y = admissions, colour = region, text = text)) +
        geom_point(alpha = 0.7, size = 1.5) +
        geom_line(aes(group = region), alpha = 0.7, size = 1.5) +
        scale_color_viridis(discrete = T, name = "Region") +
        theme_ipsum() +
        labs(x = "Date", y = "Number of Admissions")

ggpp <- ggplotly(ggp, tooltip = "text")

## if(sus_seb_combination == 0) {
##     date.adm.str <- date.adm_sus
## } else {
##     date.adm.str <- date.adm_seb
## }

plot.filename <- file.path(dirname(tmpFile), paste0("adm_plot", date.adm.str, ".jpg"))
ggsave(plot.filename,
       ggp,
       width = 9.15,
       height = 6)
saveWidget(ggpp, file=file.path(dirname(tmpFile), paste0("adm_plot", date.adm.str, ".html")))
