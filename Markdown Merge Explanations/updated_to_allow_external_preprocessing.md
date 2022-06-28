## How to run with preprocessed data

1) Checkout the changes from the branch SMED_Test_JK_Branch or Add_4_doses to the following files:

    - config.R
        - Modified settings and run names
    - run.R
        - Modified to ensure file names are correct and Angelos' pre-processing works and correct filenames are initialised
    - R/data/format_sero.R
        - Modified serology to filter to the right cutoff and reads in preprocessed serology data
    - R/data/format_deaths.R
        - Modified deaths to filter to the right cutoff and reads in preprocessed deaths data
    - set_up_inputs.R
        - Modified code to ensure the correct end dates are selected


2) To run the model the differences in the config file will be:

    -  ```r 
        ## Serology flags
        Use_preprocessed_serology <- T
        preprocessed_sero_date <- ymd("20220627") # specify the date used for the preprocessed data

        sero_cutoff_flag <- T # Modified to true to run with serology cutoff
        
        if(sero_cutoff_flag) {
            serology.delay <- 25L
            sero.end.date <- ymd(20211030)  # removed addtion of serology delay to this term
        }
        ```
        - The sero.end.date should be the last date of the data due to the later transform using the delay
    - ```r
        ## Deaths Flags
        Use_preprocessed_deaths <- T
        ```
    - Note make sure `date.data` matches the preprocessed date
