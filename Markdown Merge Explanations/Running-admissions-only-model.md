## How to run admissions only model

1) Download previous_run_input.zip from previous email
2) Extract previous_run_input.zip into {repo}/data
    
    - Note that this may overwrite the previous data in here so check that any runs stored manually may need to be backed up
    - This shouldn't be necessary as the zip file will contain the data used for a run for all-hospitalisations and one for admissions only

3) Checkout the changes from the branch SMED_Test_JK_Branch to the following files:

    - Config.R (in particular the following):
        ```r
        # Variable to determine whether or not the admissions (T) or admissions + diagnoses (F) should be used
                # Initially uses a Naive implementation
        admissions_only.flag <- T
        ## ## Value to note which combination of hospital data to use sus (0), sus + sebs (1), sebs only (2) or sus (preprocessed) + sebs (3)
        sus_seb_combination <- 3L
        ## ##Value to note how many days to remove from the end of the dataset
        adm_sus.strip_days <- 30L
        adm_seb.strip_days <- 2L
        seb_report_delay <- 1L  ## Used within this file, so can't be moved.
        date.adm_seb <- ymd(20220401)
        date.adm_sus <- ymd(20210930)
        date.adm.str <- lubridate::as_date(ifelse(sus_seb_combination > 0,
                                                          date.adm_seb - adm_seb.strip_days,
                                                          date.adm_sus - adm_sus.strip_days))

        ## ## file.locs for admissions for geography linkers (with colname links)
        adm.sus.geog_link.loc <- "utility_files/lad_to_region.csv"
        adm.sus.geog_link <- "LAD19CD"
        adm.sus.region_col <- "RGN19NM"
        adm.seb.geog_link.loc <- "utility_files/trust lookup for paul.xlsx"
        adm.seb.geog_link <- "Trust_code"
        adm.seb.region_col <- "phec_nm"

        ## ## File names of pre-processed SUS data if it is to be used.
        if(!admissions_only.flag) { ## Settings for preprocessed all hospitalisations
            preprocessed_sus_names <- paste0("2022-01-02_", regions, "_6ag_counts.txt")
            sus_old_tab_sep <- T
            preprocessed_sus_csv_name <- "admissions_data_all_hosp.csv"
        } else { ## Settings for admissions_only
            preprocessed_sus_names <- paste0("2022-03-09_", regions, "_6ag_counts.txt")
            sus_old_tab_sep <- F
            preprocessed_sus_csv_name <- "admissions_data_admissions_only.csv"
        }

        names(preprocessed_sus_names) <- regions
        print(preprocessed_sus_names)

        ## ## Admissions flags/dates
        ## adm.end.date <- date.data - adm_seb.strip_days ## Set this value if we want to truncate the data before its end.
        # adm_sus.end.date <- ymd(20210505) ## New date breakpoint
        adm_sus.end.date <- ymd(20201014)
        ```
    - format_hosp_admissions.R
    - run.R
        - Changes here only affect when code is preprocessed separately to the postprocessing
        - In particular the changes are in the following lines:
            ```r
            if(adm.flag){
                if(!format.inputs) {
                    adm_csv_fname <- ifelse(admissions_only.flag, "admissions_data_admissions_only.csv", "admissions_data_all_hosp.csv")
                    adm.sam <- read_csv(file.path(data.dirs["adm"], adm_csv_fname))
                    file.copy(file.path(data.dirs["adm"], adm_csv_fname), out.dir)
                    file.rename(file.path(data.dirs["adm"], adm_csv_fname), file.path(data.dirs["adm"], "admissions_data.csv"))
                    admsam.files <- paste0(data.dirs["adm"], "/", date.adm.str, "_", regions, "_", nA_adm, "ag_counts.txt")
                }   
            }
            ```

4) To run the model using the admissions only model the differences in the config file will be:
    -  ```r 
        admissions_only_flag <- T
        ```
    - ```r
        sus_seb_combination <- 3L ## Note that this can be used for both (i.e. all hospitalisations as well, but will use a previous run instead of directly reading in SUS data)
        ```