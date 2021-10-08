suppressMessages(require(lubridate))
suppressMessages(require(tidyverse))
suppressMessages(require(gsubfn))
suppressMessages(require(readxl))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################
## Where to find the data, if NULL use command line argument
if(!exists("sero.loc")){ ## Set to default format for the filename
    #input.loc <- build.data.filepath(subdir = "raw", "serology")
    input.loc <- "data/raw/serology"
    input.loc <- "~/CoVID-19/Data streams/Serology"
    ## List the possible files in the directory
    sero.loc <- file.info(file.path(input.loc, list.files(path=input.loc, pattern=glob2rx("*xlsx|*csv"))))
    ## Pick the most recently added
    input.loc <- rownames(sero.loc)[which.max(sero.loc$mtime)]
} else {
    if(is.null(sero.loc)){
        input.loc <- commandArgs(trailingOnly = TRUE)[1]
    } else {
        if(startsWith(sero.loc, "/")) input.loc <- sero.loc
        else input.loc <- build.data.filepath(subdir = "raw", "serology", sero.loc)
    }
}
csv.flag <- grepl("*csv$", input.loc)

## Set the max and min dates for the data
earliest.date <- start.date
latest.date <- sero.end.date

## What is the date of publication of these data? If not specified, try to extract from filename
if(!exists("date.sero")){
    fl.name <- basename(input.loc)
    date.sero.str <- strapplyc(fl.name, "[0-9]{8,}", simplify = TRUE)
    if(length(date.sero.str[[1]]) == 0) date.sero.str <- sero.end.date
}

## Define an age-grouping
if(!exists("age.agg")){
    age.agg <- c(0, 1, 5, 15, 25, 45, 65, 75, Inf)
    age.labs <- c("<1yr","1-4","5-14","15-24","25-44","45-64","65-74", "75+")
    }
nA <- length(age.labs)

if(!exists("regions")){
    regions <- c("London", "Outside_London")
}

## Map our names for columns (LHS) to data column names (RHS)
possible.col.names <- list(
    surv = c("surv", "Collection"),
    age = "age",
    region = c("region", "Region", "NHS_Region"),
    sample_date = c("sampledate", "SampleDate"),
    Eoutcome = c("EuroImm_outcome", "EuroImmun_outcome", "euroimmun_outcome"),
    Eresult = c("EuroImmun_units", "EuroImm_Units", "euroimmun_units"),
    Routcome = c("RBD_outcome", "RBD_outcome", "rbd_outcome"),
    Rresult = c("RBD_units", "RBD_Units", "rbd_units"),
    RNoutcome = c("RocheN_outcome", "roche_n_outcome", "Roche_N_outcome"),
    RNresult = c("RocheN_units", "roche_n_units", "Roche_N_units"),
    RSoutcome = c("RocheS_outcome", "roche_s_outcome", "Roche_S_outcome"),
    RSresult = c("RocheS_units", "roche_s_units", "Roche_S_units")
)
if(region.type == "ONS") possible.col.names$ONS_region = "ONS_Region"

## names(read_xlsx(input.loc,sheet = 1, n_max = 0))

if(csv.flag){
    input.col.names <- suppressMessages(names(read_csv(input.loc, n_max = 0)))
} else input.col.names <- suppressMessages(names(read_xlsx(input.loc,sheet = 1, n_max = 0)))
is.valid.col.name <- function(name) {name %in% input.col.names}
first.valid.col.name <- function(names) {first.where.true(names, is.valid.col.name)}
col.names <- lapply(possible.col.names, first.valid.col.name)
invalid.col.names <- sapply(col.names, is.null)
if (any(invalid.col.names)) {
	names.invalid.cols <- paste0(names(possible.col.names)[invalid.col.names], collapse = ", ")
	stop(paste("No valid column name for:", names.invalid.cols))
}

## Given a row in the sero data file, return its region, formatted with no spaces
if (region.type == "NHS") {
	get.region <- function(x) str_replace_all(x$region, " ", "_")
} else {
	get.region <- function(x) str_replace_all(x$ONS_region, " ", "_")
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

if(!exists("serosam.files")){
    serosam.files <- build.data.filepath("RTM_format/serology",
                                         "sero",
                                         date.sero.str,
                                         "_",
                                         regions,
                                         "_",
                                         nA,
                                         "ages_samples.txt")
    seropos.files <- build.data.filepath("RTM_format/serology",
                                         "sero",
                                         date.sero.str,
                                         "_",
                                         regions,
                                         "_",
                                         nA,
                                         "ages_positives.txt")
}
## Which columns are we interested in?
sero.col.args <- list()
sero.col.args[[col.names[["surv"]]]] <- col_character()
sero.col.args[[col.names[["age"]]]] <- col_integer()
sero.col.args[[col.names[["region"]]]] <- col_character()
sero.col.args[[col.names[["sample_date"]]]] <- col_date(sero.date.fmt)
sero.col.args[[col.names[["Eoutcome"]]]] <- col_character()
sero.col.args[[col.names[["Eresult"]]]] <- col_double()
sero.col.args[[col.names[["Routcome"]]]] <- col_character()
sero.col.args[[col.names[["Rresult"]]]] <- col_double()
sero.col.args[[col.names[["RNoutcome"]]]] <- col_character()
sero.col.args[[col.names[["RNresult"]]]] <- col_double()
sero.col.args[[col.names[["RSoutcome"]]]] <- col_character()
sero.col.args[[col.names[["RSresult"]]]] <- col_double()
if(region.type == "ONS") sero.col.args[[col.names[["ONS_region"]]]] <- col_character()
sero.cols <- do.call(cols_only, sero.col.args)

get.col.type <- function(x){
    ifelse(identical(x, col_character()), "text",
    ifelse(identical(x, col_double()), "numeric",
    ifelse(identical(x, col_integer()), "numeric",
    ifelse(identical(x, col_date(sero.date.fmt)), "date", "skip"))))
}

## Reading in the data ##
print(paste("Reading from", input.loc))
strPos <- c("+", "Positive", "positive")
#! read xlsx line needs to have improved robustness


if(csv.flag){
    sero.dat <- read_csv(input.loc, col_types = sero.cols)
} else {
    fields <- sapply(sero.col.args[input.col.names], get.col.type)
    sero.dat <- read_xlsx(input.loc, sheet = 1, col_types = fields)
}
sero.dat <- sero.dat %>% 
    rename(!!!col.names) %>%
    pivot_longer(contains("outcome"), names_to = "assay", values_to = "outcome") %>% 
    filter(!is.na(outcome)) %>% 
    mutate(SDate = fuzzy_date_parse(sample_date),
           Positive = endsWith(outcome, "+") | outcome == "Positive")

## Set which collection to examine
if(NHSBT.flag) {
  sero.dat <- sero.dat %>%
    filter(startsWith(surv, "NHSBT") | startsWith(surv, "NHS BT"))
} else {
  sero.dat <- sero.dat %>%
      filter(startsWith(surv, "RCGP") | startsWith(surv, "rcgp") |
             ((startsWith(surv, "NHSBT") | startsWith(surv, "NHS BT")) & SDate < ymd("2020-05-22"))
             )
}

## Apply filters to get only the data we want.
sero.dat <- sero.dat %>%
    filter(!is.na(region), SDate <= sero.end.date) %>%
    mutate(region = get.region(.),
           age.grp = cut(age, age.agg, age.labs, right = FALSE, ordered_result = T),
           date = SDate - serology.delay)

## Set correct name for the assay we want to examine
if(RocheS.flag) {
  roche <- "RSoutcome"
} else {
  roche <- "RNoutcome"
}

## ## Get into format for use by the rtm
## First denominators
rtm.sam <- sero.dat %>%
    filter((assay == "Eoutcome" & as_date(date) < ymd("2020-05-22")) | (assay == roche & as_date(date) >= ymd("2020-05-22"))) %>% 
    group_by(date, region, age.grp, .drop = FALSE) %>%
    tally %>%
    right_join(expand_grid(date = as_date(earliest.date:latest.date), ## Add unsamples region, date and age group combinations
                           region = regions,
                           age.grp = age.labs),
               by = c("date", "region", "age.grp")
               ) %>%
    replace_na(list(n = 0)) %>% ## 0 if just added
    arrange(date)
## Then positives
rtm.pos <- sero.dat %>%
    filter((assay == "Eoutcome" & as_date(date) < ymd("2020-05-22")) | (assay == roche & as_date(date) >= ymd("2020-05-22"))) %>% 
    filter(Positive) %>%
    group_by(date, region, age.grp, .drop = FALSE) %>%
    tally %>%
    right_join(expand_grid(date = as_date(earliest.date:latest.date), ## Add unsamples region, date and age group combinations
                           region = regions,
                           age.grp = age.labs),
               by = c("date", "region", "age.grp")
               ) %>%
    replace_na(list(n = 0)) %>% ## 0 if just added
    arrange(date)

## Write rtm data outputs to file
#serosam.files <- str_replace_all(serosam.files, date.data, date.sero.str)
#seropos.files <- str_replace_all(seropos.files, date.data, date.sero.str)
names(serosam.files) <- names(seropos.files) <- regions

# If directory required doesn't exist create it
walk(dirname(serosam.files), ~dir.create(., showWarnings = F, recursive = T))

for(reg in regions){
    region.sam <- pivot_wider(rtm.sam %>%
                              filter(region == reg),
                              id_cols = 1,
                              names_from = age.grp,
                              values_from = n)
    region.pos <- pivot_wider(rtm.pos %>%
                              filter(region == reg),
                              id_cols = 1,
                              names_from = age.grp,
                              values_from = n)
    
    tmpFile <- serosam.files[reg]

    print(paste("Writing to",
                tmpFile,
                "(",
                sum(region.sam[, -1]),
                "total samples,",
                nrow(region.sam),
                "rows.)"
                )
          )

    region.sam %>%
        write_tsv(
            tmpFile,
            col_names = FALSE
        )

    tmpFile <- seropos.files[reg]

    print(paste("Writing to",
                tmpFile,
                "(",
                sum(region.pos[, -1]),
                "total positives",
                nrow(region.pos),
                "rows.)"
                )
          )
    
    region.pos %>%
        write_tsv(
            tmpFile,
            col_names = FALSE
            )

    }

# If directory required doesn't exist create it
if(!file.exists(out.dir))
    dir.create(out.dir, recursive = T)

## Save the data
write_csv(rtm.sam, file.path(out.dir, "sero_samples_data.csv"))
write_csv(rtm.pos, file.path(out.dir, "sero_positives_data.csv"))

## Save a quick plot of the data..
require(ggplot2)
require(viridis)
require(hrbrthemes)
require(plotly)
require(htmlwidgets)

strAssay <- ifelse(RocheS.flag, "S", "N")
strCollection <- ifelse(NHSBT.flag, "NHSBT", "RCGP")

rtm.plot <- rtm.sam %>%
    inner_join(rtm.pos %>% rename(p = n)) %>%
    group_by(date, region) %>%
    summarise(N = sum(n),
              P = sum(p)) %>%
    mutate(X = P / N) %>%
    filter(!is.nan(X)) %>%
    mutate(text = paste(region, ": N = ", N))

gp <- ggplot(rtm.plot, aes(x = date, y = X, size = N, colour = region, text = text)) +
    geom_point(alpha = 0.7, position = position_dodge2(width = 0.6, preserve = "single")) +
    scale_size(range = c(1, 10), name="Sample size, N") +
    scale_color_viridis(discrete = TRUE, name = "Region") +
    theme_ipsum() +
    ylab("Proportion positive")

pp <- ggplotly(gp, tooltip = "text")

plot.filename <- file.path(dirname(tmpFile), paste0("sero_plot", date.sero.str, ".jpg"))
ggsave(plot.filename,
       gp,
       width = 9.15,
       height = 6,
       title = glue::glue("Assay: Roche-", strAssay, "; Collection: ", strCollection))
saveWidget(pp, file=file.path(dirname(tmpFile), paste0("sero_plot", date.sero.str, strCollection, "_", strAssay, ".html")))

## pp.focus <- ggplotly(gp+ylim(c(0, 0.25)), tooltip = "text")

## plot.filename <- file.path(dirname(tmpFile), paste0("sero_plot_focus", date.sero.str, ".jpg"))
## ggsave(plot.filename,
##        gp,
##        width = 9.15,
##        height = 6,
##        main = glue::glue("Assay: Roche-", ifelse(RocheS.flag, "S", "N")"; Collection: ", ifelse(NHSBT.flag, "NHSBT", "RCGP")))
## saveWidget(pp.focus, file=file.path(dirname(tmpFile), paste0("sero_plot_focus", date.sero.str, ".html")))

rtm.plot.by.age <- rtm.sam %>%
    inner_join(rtm.pos %>% rename(p = n)) %>%
    mutate(X = p / n) %>%
    filter(!is.nan(X)) %>%
    mutate(text = paste(region, ": N = ", n))

gap <- ggplot(rtm.plot.by.age, aes(x = date, y = X, size = n, colour = region, text = text)) +
    geom_point(alpha = 0.7, position = position_dodge2(width = 0.6, preserve = "single")) +
    scale_size(range = c(1, 10), name = "Sample size, N") +
    scale_colour_viridis(discrete = TRUE, name = "Region") +
    theme_ipsum() +
    ylab("Proportion positive") +
    facet_wrap(~age.grp)
pap <- ggplotly(gap, tooltip = "text")

plot.filename <- file.path(dirname(tmpFile), paste0("sero_plot_byage", date.sero.str, ".jpg"))
ggsave(plot.filename,
       gap,
       width = 9.15,
       height = 6,
       title = glue::glue("Assay: Roche-", strAssay, "; Collection: ", strCollection))
saveWidget(pap, file=file.path(dirname(tmpFile), paste0("sero_plot", date.sero.str, strCollection, "_", strAssay, ".html")))

stop()
