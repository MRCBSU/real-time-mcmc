## Where to find the data as downloaded from REACT?
## The data come from .csv files corresponding to each round of the study
## separately for the NHS and the GO regions.
## The data on the rounds have also arrived in batches, with differing levels of temporal aggregation. May use earlier versions of the data for finer temporal separation.

if(!exists("react.loc")) ## Choose the default location for these data
    react.loc <- build.data.filepath(subdir = "raw", "react")

## Number of rounds of the survey so far should be defined in config.R
## Use this and the geography to piece together the data file names.
react.loc <- file.path(react.loc, paste("round", 1:react.rounds, paste0(tolower(region.type), ".csv"), sep = "_"))

## Map our names for columns (LHS) to data column names (RHS)
possible.col.names <- list(
    age = "age_group",
    region = c("nhs_region", "region"),
    sdate = c("date", "react_week_start_date"),
    npos = "number_positive",
    nsam = "number_samples"
)

get.age.groupings <- function(x){
    
    }

## Set output object
react.dat <- tibble()

## Loop through the rounds appending each dataset to the above output object
for(i in 1:react.rounds)
{

    ## Check the dataset to be read in has the required columns
    input.col.names <- suppressMessages(names(read_csv(react.loc[i], n_max = 0)))
    is.valid.col.name <- function(name) {name %in% input.col.names}
    first.valid.col.name <- function(names) {first.where.true(names, is.valid.col.name)}
    col.names <- lapply(possible.col.names, first.valid.col.name)
    invalid.col.names <- sapply(col.names, is.null)
    if(any(invalid.col.names)){
        names.invalid.cols <- paste0(names(possible.col.names)[invalid.col.names], collapse = ", ")
        stop(paste("No valid column name for:", names.invalid.cols, "in REACT round", i))
    }

    ## How are the columns defined in the data we will read in?
    react.col.args <- list()
    react.col.args[[col.names[["age"]]]] <- col_character()
    react.col.args[[col.names[["region"]]]] <- col_character()
    react.col.args[[col.names[["sdate"]]]] <- col_date()
    react.col.args[[col.names[["npos"]]]] <- col_integer()
    react.col.args[[col.names[["nsam"]]]] <- col_integer()
    react.cols <- do.call(cols, react.col.args)

    ## Reading in the data ##
    print(paste("Reading from", react.loc[i]))

    react.dat <- bind_rows(react.dat,
                           read_csv(react.loc[i],
                                    col_types = react.cols) %>%
                           rename(!!!col.names) %>%
                           select(!!!(names(col.names))) %>%
                           mutate(sdate = fuzzy_date_parse(sdate),
                                  round = i) %>%
                           )
}

## ## Outstanding issues - 
## remove NA regions and rename appropriately
react.dat <- react.dat %>%
    filter(!is.na(region)) %>%
    get.region()
## add on half a week for the later rounds.
## aggregate to the correct age groups


## Need to test this with the ONS.
