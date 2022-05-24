## How to use the new run interface

1) Checkout the folder Utils from the root of the directory from the branch SMED_Test_JK_Branch or Add_4_doses
2) From any folder you would like to run one of scripts (i.e. from the root for the preprocessing or from the  model_runs/date for the main running or postprocessing)
3) Follow the instructions in the shell to run the appropriate script
    - Note: that this will let you customise the array inputs where appropriate
    - Note: Other parameters will need ot be edited in the appropriate scripts such as config.R or the submission script as these are called by the interface

## How to run cutoff models

1) Checkout the changes from the branch SMED_Test_JK_Branch or Add_4_doses to the following files:

    - config.R
        - Modified settings and run names
    - run.R
        - Modified to ensure file names are correct and Angelos' pre-processing works
    - R/data/format_sero.R
        - Modified serology to filter to the right cutoff and write to custom filenames
    - R/data/format_deaths.R
        - Modified deaths to filter to the right cutoff and write to custom filenames
    - R/data/format_hosp_admissions.R
        - Modified deaths to filter to the right cutoff and write to custom filenames
    - set_up_inputs.R
        - Modified code to ensure the correct end dates are selected


2) To run the model using the admissions only model the differences in the config file will be:

    -  ```r 
        ## Serology flags
        sero_cutoff_flag <- T # Modified to true to run with serology cutoff
        
        if(sero_cutoff_flag) {
            serology.delay <- 25L
            sero.end.date <- ymd(20211030) + serology.delay # Modified to choose cutoff date
        }
        ```
        - The sero.end.date should be the last date of the data added to the serology delay due to the later transform using the delay
    - ```r
        ## Deaths Flags
        use_deaths_up_to_now_flag <- F ## modified to false to use deaths cutoff
        custom_deaths_end_date <- lubridate::ymd("20211030") + 25 ## Modified to use cutoff date
        ```
        - Note that the cutoff is applied when use_deaths_up_to_now_flag is set to False
    - ```r
        ## Admissions Flags
        # Flag to determine whether to cutoff the hospitalisation datastream early (T => use cutoff)
cutoff_hosps_early <- T
        date_early_cutoff_hosps <- ymd(20211030) ## Modified to use cutoff date
        ```
        - Note that the cutoff is applied when use_deaths_up_to_now_flag is set to False

## Addition of 4th doses to the preprocessing

1) Checkout the following files from the branch Add_4_doses
    - config.R
        - modifies the number of doses used in the model and the run name
    - run.R
        - Ensures the file names read in are correct and fixes the doses for angelos' pipeline
    - R/data/format_vaccinations.R
        - Ensure all the preprocessed files are present and backup for naming the vaccination files
    - mod_inputs.Rmd
        - Passes the filenames to the c++ code and the flags to indicate what combination of vaccine doses should be used in the model
    - mod_pars.Rmd
        - Passes the priors on the efficacy for 4 doses to the model
    - set_up_pars.R
        - Where the preprocessing determines the efficacy of the vaccine data and ensures the priors are selected for each dose

2) The latest vaccination data should be unzipped into RTM_format (with the file name something similar to 20220413_four_doses.zip)

3) In config.r make sure the following variable is set `vacc.n_doses <- 4L`
