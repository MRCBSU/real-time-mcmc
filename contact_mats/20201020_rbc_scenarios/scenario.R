library(tidyverse)
library(contactsr)

future_scenarios <- function(data, region) {
  readRDS("timeuse.rds") %>% dplyr::filter(NHSER19NM == region) -> tu_df

  startDate <- lubridate::now() %>% lubridate::as_date()

  tu_df %>% dplyr::filter(Week == lubridate::isoweek(startDate)) -> tu_df
  readRDS("current_activity_levels.rds") %>% dplyr::filter(NHSER19NM == region) -> df
  #startDate <- df$Date %>% unique()
  df %>%
    dplyr::select(-Date, -Week) -> activity_df

  base_matrices <- readRDS("split_polymod_matrices.rds")[[region]]

  # Pick all contacts (not just physical)
  base_matrices %>% purrr::map(function(lst) lst$all) -> base_matrices

  # Combine them for later
  base_matrices %>% purrr::reduce(function(acc, cm_lst) {
    acc$matrix <- acc$matrix + cm_lst$matrix
    acc
  }) -> combined_matrix

  # Create the future activity levels
  activity_df %>% 
    dplyr::inner_join(tu_df %>% 
                     dplyr::select(location, activity) %>% dplyr::distinct()) -> act_loc_df

  split(data, data$Date) %>% purrr::accumulate(.init = act_loc_df, function(acc, df) {
    loc_df <- df %>% dplyr::filter(is.na(Activity)) %>%
      dplyr::select(Location, value) %>%
      dplyr::rename(location = Location, nlvalue = value)
    act_df <- df %>% dplyr::filter(is.na(Location)) %>%
      dplyr::select(Activity, value) %>%
      dplyr::rename(activity = Activity, navalue = value)

    acc %>% dplyr::left_join(loc_df) %>%
      dplyr::left_join(act_df) %>%
      dplyr::mutate(value = dplyr::case_when(!is.na(navalue) ~ navalue,
                                             !is.na(nlvalue) ~ nlvalue,
                                             T ~ value)) %>%
      dplyr::select(-nlvalue, -navalue)
  }) -> activity_lst

  activity_lst[[as.character(startDate)]] <- activity_lst[[".init"]]

  # Get rid of ".init" and sort by name/date
  nms <- names(activity_lst)
  names(nms) <- nms
  nms %>% sort() %>%
    purrr::keep(function(nm) nm != ".init") %>%
    purrr::map(function(nm) activity_lst[[nm]]) -> activity_lst

  activity_lst %>% purrr::imap(function(act_df, date) {
    mask_lst <- contactsr::timeuse_location_matrix(tu_df, act_df %>% dplyr::select(-location))

    if (!all(sort(names(mask_lst)) == sort(names(base_matrices))))
      stop("Error detected in locations")

    names(mask_lst) %>% purrr::map(function(name) {
      base_matrices[[name]]$matrix*mask_lst[[name]]
    }) %>% purrr::reduce(function(acc, m) acc + m) -> mat

    mat/combined_matrix$matrix
  })
}

# PAUL: Below here is how you define the scenarios and create the (relative) matrices over time

# You need to have an up to date version of contactsr installed
#
# devtools::install_git("https://gitlab.com/epidemics-r/contactsr")

# Define your scenario by the Date the changes happen, 
# which location changes (home, school, work, leisure, transport, otherplace)
# or which activity (school, university, parks, etc. (see below for a complete list)).
# and what the new value is (0-1)
#
# DO N0T define both Location and Activity (one of them should always be NA)
# For your use case it might be easier to read this in from a csv file, than to define the table in the code
# changes_df <- dplyr::tibble(Date = c("2020-10-23", "2020-10-23", "2020-11-29", "2020-11-29"), 
#                            Location = c(NA, "leisure", "leisure", NA),
#                            Activity = c("school", NA, NA, "university"),
#                            value = c(0.85, 0.5, 1, 1))

# Activities in general are more fine grained than locations
# For example the school location is split into activities school (primary/secondary school) and university
# Complete list of possible locations and activities 
#readRDS("timeuse.rds") %>% dplyr::select(location, activity) %>%
#  dplyr::distinct()

#future_scenarios(changes_df, "England")

# Should be made into a function
simple_scalar <- function(scalar) {
  readRDS("timeuse.rds") %>% dplyr::select(location, activity) %>%
    dplyr::distinct() -> df

  all_activities <- df$activity %>% unique()
  acs <- setdiff(all_activities, c("home", "school"))

  dplyr::tibble(Activity = acs) %>%
    dplyr::mutate(Date = "2020-10-25",
                  Location = NA_character_, value = scalar) -> changes_df
  future_scenarios(changes_df, "England") -> lst

  message("Relative dominant eigenvalue ", 
          1/(eigen(lst[[1]], only.values = T)$value[[1]]/eigen(lst[[2]], only.values = T)$value[[1]]))

  lst[[2]]
}

simple_scalar(0) -> value_00
## write.csv(value_00, row.names = F, file = "home_school_00.csv")
simple_scalar(0.2) -> value_02
## write.csv(value_02, row.names = F, file = "home_school_20.csv")
