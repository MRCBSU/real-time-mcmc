suppressMessages(library(lubridate))
suppressMessages(library(tidyverse))
suppressMessages(library(gsubfn))

#########################################################
## Inputs that should (or may) change on a daily basis
#########################################################
## Where to find the data, if NULL use command line argument
if(!exists("sero.loc")){ ## Set to default format for the filename
    #input.loc <- build.data.filepath(subdir = "raw", "serology")
    input.loc <- "data/raw/serology"
    ## List the possible files in the directory
    sero.loc <- file.info(file.path(input.loc, list.files(path=input.loc, pattern=glob2rx("202*csv"))))
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

earliest.date <- start.date
latest.date <- sero.end.date

## What is the date of publication of these data? If not specified, try to extract from filename
if(!exists("date.sero")){
    fl.name <- basename(input.loc)
    date.sero.str <- strapplyc(fl.name, "[0-9]{8,}", simplify = TRUE)
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
    surv = "surv",
    age = "age",
    region = c("region", "Region"),
    sample_date = c("sampledate", "SampleDate"),
    Eoutcome = c("EuroImm_outcome", "EuroImmun_outcome", "euroimmun_outcome"),
    Eresult = c("EuroImmun_units", "EuroImm_Units", "euroimmun_units"),
    Routcome = c("RBD_units", "RBD_outcome", "rbd_outcome"),
    Rresult = c("RBD_units", "RBD_Units", "rbd_units")
)
if(region.type == "ONS") possible.col.names$ONS_region = "ONS_Region"

input.col.names <- suppressMessages(names(read_csv(input.loc, n_max=0)))
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
sero.col.args[[col.names[["sample_date"]]]] <- col_character()
sero.col.args[[col.names[["Eoutcome"]]]] <- col_character()
sero.col.args[[col.names[["Eresult"]]]] <- col_double()
sero.col.args[[col.names[["Routcome"]]]] <- col_character()
sero.col.args[[col.names[["Rresult"]]]] <- col_double()
if(region.type == "ONS") sero.col.args[[col.names[["ONS_region"]]]] <- col_character()
sero.cols <- do.call(cols_only, sero.col.args)

## Reading in the data ##
print(paste("Reading from", input.loc))
strPos <- c("+", "Positive", "positive")
sero.dat <- read_csv(input.loc,
                     col_types = sero.cols) %>%
    rename(!!!col.names) %>%
    mutate(SDate = fuzzy_date_parse(sample_date),
           Positive = endsWith(Eoutcome, "+") | Eoutcome == "Positive")

## Apply filters to get only the data we want.
sero.dat <- sero.dat %>%
    filter(startsWith(surv, "NHSBT")) %>%
    filter(!is.na(region), SDate <= sero.end.date) %>%
    mutate(region = get.region(.),
           age.grp = cut(age, age.agg, age.labs, right = FALSE, ordered_result = T),
           date = SDate - serology.delay)

## ## Get into format for use by the rtm
## First denominators
rtm.sam <- sero.dat %>%
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

## Save the data
write_csv(rtm.sam, file.path(out.dir, "sero_samples_data.csv"))
write_csv(rtm.pos, file.path(out.dir, "sero_positives_data.csv"))

## Save a quick plot of the data..
require(ggplot2)
require(viridis)
require(hrbrthemes)
require(plotly)
require(htmlwidgets)

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
       height = 6)
saveWidget(pp, file=file.path(dirname(tmpFile), paste0("sero_plot", date.sero.str, ".html")))

pp.focus <- ggplotly(gp+ylim(c(0, 0.25)), tooltip = "text")

plot.filename <- file.path(dirname(tmpFile), paste0("sero_plot_focus", date.sero.str, ".jpg"))
ggsave(plot.filename,
       gp,
       width = 9.15,
       height = 6)
saveWidget(pp.focus, file=file.path(dirname(tmpFile), paste0("sero_plot_focus", date.sero.str, ".html")))
