## How to run cutoff models

1) Checkout the changes from the branch SMED_Test_JK_Branch to the following files:

    - config.R
        - Modified settings and run names
    - run.R
        - Modified to ensure file names are correct and Angelos' pre-processing works
    - R/data/format_sero.R
        - Modified serology to filter to the right cutoff and write to custom filenames
    - R/data/format_deaths.R
        - Modified deaths to filter to the right cutoff and write to custom filenames


2) To run the model using the admissions only model the differences in the config file will be:

    -  ```r 
        sero_cutoff_flag <- T # Modified to true to run with serology cutoff
        if(sero_cutoff_flag) {
            serology.delay <- 25
            sero.end.date <- ymd(20211030) # Modified to choose cutoff date
        }
        ```
    - ```r
        ## Deaths Flags
        use_deaths_up_to_now_flag <- F ## modified to false to use deaths cutoff
        custom_deaths_end_date <- lubridate::ymd("20211030") ## Modified to use cutoff date
        ```
