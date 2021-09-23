suppressMessages(require(lubridate))
suppressMessages(require(tidyverse))
suppressMessages(require(readxl))
suppressMessages(require(gsubfn))

#########################################################
## Inputs that should (or may) change on a weekly basis
#########################################################
## Where to find the data, if NULL use command line argument
if(!exists("adm.loc")){ ## Set to default format for the filename
    #input.loc <- build.data.filepath(subdir = "raw", "serology")
    input.loc <- "data/raw/admissions"
    ## List the possible files in the directory
    adm.loc <- file.info(file.path(input.loc, list.files(path=input.loc, pattern=glob2rx("202*xlsx"))))
    ## Pick the most recently added
    input.loc <- rownames(adm.loc)[which.max(adm.loc$mtime)]
} else {
    if(is.null(adm.loc)){
        input.loc <- commandArgs(trailingOnly = TRUE)[1]
    } else {
        if(startsWith(adm.loc, "/")) input.loc <- adm.loc
        else input.loc <- build.data.filepath(subdir = "raw", "admissions", sero.loc)
    }
}

## What is the date of publication of these data? If not specified, try to extract from filename
if(!exists("date.adm")){
    fl.name <- basename(input.loc)
    date.adm.str <- strapplyc(fl.name, "[0-9]{8,}", simplify = TRUE)
    date.adm <- as_date(paste(substr(date.adm.str,1,4),
                                     substr(date.adm.str,5,6),
                                     substr(date.adm.str,7,8),
                                     sep = "-"))
}

earliest.date <- start.date
latest.date <- date.adm - 1

## Define an age-grouping
if(!exists("age_adm.agg")){
    age_adm.agg <- c(0, 6, 18, 65, 85, Inf)
    age_adm.oldlabs <- c("0_5", "6_17", "18_64", "65_84", "85", "")
    age_adm.labs <- c("0-5", "6-17", "18-64", "65-84", "85+", "NA") ## "All ages"
}
nA_adm <- length(age_adm.labs)

#! This behavioiur is not yet functional
# if(!exists("regions")){
#     regions <- c("London", "Outside_London")
# }

possible.col.names <- list(
    region = c("region", "Region", "Geography", "NHS_Region"),
    date = c("DateVal"),
    report_level = c("ReportLevel"),
    trust_code = c("TrustCode"),
    trust_name = c("TrustName"),
    admissions = c("hospital_inc"),
    admissions0_5 = c("hospital_inc_0-5"),
    admissions6_17 = c("hospital_inc_6-17"),
    admissions18_64 = c("hospital_inc_18-64"),
    admissions65_84 = c("hospital_inc_65-84"),
    admissions85 = c("hospital_inc_85+"),
    diagnoses = c("hospital_inc_new"),
    diagnoses0_5 = c("hospital_inc_new_0-5"),
    diagnoses6_17 = c("hospital_inc_new_6-17"),
    diagnoses18_64 = c("hospital_inc_new_18-64"),
    diagnoses65_84 = c("hospital_inc_new_65-84"),
    diagnoses85 = c("hospital_inc_new_85+")
)

input.col.names <- suppressMessages(names(read_xlsx(input.loc, sheet = "Extracted Data", n_max =0)))
is.valid.col.name <- function(name) {name %in% input.col.names}
first.valid.col.name <- function(names) {first.where.true(names, is.valid.col.name)}
col.names <- lapply(possible.col.names, first.valid.col.name)
invalid.col.names <- sapply(col.names, is.null)
if (any(invalid.col.names)) {
	names.invalid.cols <- paste0(names(possible.col.names)[invalid.col.names], collapse = ", ")
	stop(paste("No valid column name for:", names.invalid.cols))
}

if (region.type == "NHS") {
	get.region <- function(x) str_replace_all(x$region, " ", "_")
} else {
	get.region <- function(x) str_replace_all(x$ONS_region, " ", "_")
    stop("ONS data not yet implemented")
}

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
    admsam.files <- build.data.filepath("RTM_format/admission",
                                         "adm",
                                         date.adm.str,
                                         "_",
                                         regions,
                                         "_",
                                         nA_adm,
                                         "ages_samples.txt")
    admpos.files <- build.data.filepath("RTM_format/admission",
                                         "adm",
                                         date.adm.str,
                                         "_",
                                         regions,
                                         "_",
                                         nA_adm,
                                         "ages_positives.txt")
}


adm.date.fmt <- "%y-%b-%d" #! Not right but irrelevant
adm.col.args <- list()
adm.col.args[[col.names[["region"]]]] <- col_character()
adm.col.args[[col.names[["date"]]]] <- col_date(adm.date.fmt)
adm.col.args[[col.names[["report_level"]]]] <- col_character()
adm.col.args[[col.names[["trust_code"]]]] <- col_character()
adm.col.args[[col.names[["trust_name"]]]] <- col_character()
adm.col.args[[col.names[["admissions"]]]] <- col_integer()
adm.col.args[[col.names[["admissions0_5"]]]] <- col_integer()
adm.col.args[[col.names[["admissions6_17"]]]] <- col_integer()
adm.col.args[[col.names[["admissions18_64"]]]] <- col_integer()
adm.col.args[[col.names[["admissions65_84"]]]] <- col_integer()
adm.col.args[[col.names[["admissions85"]]]] <- col_integer()
adm.col.args[[col.names[["diagnoses"]]]] <- col_integer()
adm.col.args[[col.names[["diagnoses0_5"]]]] <- col_integer()
adm.col.args[[col.names[["diagnoses6_17"]]]] <- col_integer()
adm.col.args[[col.names[["diagnoses18_64"]]]] <- col_integer()
adm.col.args[[col.names[["diagnoses65_84"]]]] <- col_integer()
adm.col.args[[col.names[["diagnoses85"]]]] <- col_integer()
adm.cols <- do.call(cols_only, adm.col.args)

get_col_type <- function(x) {
    ifelse(
        identical(x, col_character()), "text",
        ifelse(identical(x, col_double()), "numeric",
        ifelse(identical(x, col_integer()), "numeric",
        ifelse(identical(x, col_date(adm.date.fmt)), "text",
        "skip"
    ))))
}
fields <- sapply(adm.col.args[input.col.names], get_col_type)

print(paste0("Reading in data from ", input.loc))

adm.dat <- read_xlsx(input.loc, sheet = "Extracted Data", col_types = fields) %>%
    rename(!!!col.names)  %>%
    pivot_longer(where(is.numeric), names_to = "ages", values_to = "count") %>%
    filter(!is.na(count)) %>%
    mutate(type = str_extract(ages, "^(admissions|diagnoses)"),
           ages = str_remove(ages, "^(admissions|diagnoses)"),
           date = ymd(date))

adm.dat <- adm.dat %>%
    mutate(region = get.region(.),
    ages = factor(ages, levels = age_adm.oldlabs))

levels(adm.dat$ages) <- age_adm.labs

if(region.type == "NHS") {
    adm.sam <- adm.dat %>%
        filter(report_level == "Region") %>%
        select(-report_level, -trust_name, -trust_code) %>%
        right_join(crossing(date = as_date(earliest.date:latest.date),
                            region = regions,
                            ages = age_adm.labs,
                            type = c("admissions", "diagnoses")), by = c("date" = "date", "region" = "region", "ages" = "ages", "type" = "type")) %>%
        replace_na(list(count = 0)) %>%
        group_by(date, region, ages) %>%
        summarise(across(count, .fn = sum))
} else {
    #! ONS regions needs to be sorted
    adm.sam <- tribble(~data, NA)
}

## Write rtm data outputs to file
names(admsam.files) <- names(admpos.files) <- regions

## Create missing directories
walk(dirname(admsam.files), ~dir.create(., showWarnings = F, recursive = T))

# Loop over regions and save the data for them
for(reg in regions) {
    region.sam <- pivot_wider(adm.sam %>%
                              filter(region == reg),
                              id_cols = c(date, region),
                              names_from = ages,
                              values_from = count
                              )

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

# Save the data
write_csv(adm.sam, file.path(out.dir, "admissions_data.csv"))


## Create a quick plot of the data
require(ggplot2)
require(viridis)
require(hrbrthemes)
require(plotly)
require(htmlwidgets)

adm.plot <- adm.sam %>%
    group_by(date, region) %>%
    summarise(across(count, .fn = sum)) %>%
    mutate(text = paste("Region:", region, "\nAdmssions:", count, "\nDate:", date))

ggp <- ggplot(adm.plot, aes(x = date, y = count, colour = region, text = text)) +
        geom_point(alpha = 0.7, size = 1.5) +
        scale_color_viridis(discrete = T, name = "Region") +
        theme_ipsum() +
        labs(x = "Date", y = "Number of Admissions")

ggpp <- ggplotly(ggp, tooltip = "text")

plot.filename <- file.path(dirname(tmpFile), paste0("adm_plot", date.adm.str, ".jpg"))
ggsave(plot.filename,
       ggp,
       width = 9.15,
       height = 6)
saveWidget(ggpp, file=file.path(dirname(tmpFile), paste0("adm_plot", date.adm.str, ".html")))