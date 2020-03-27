# real-time-mcmc

Hosting of the PHE real-time model for pandemic influenza

## Process overview

Running the model requires the following steps, each explained in more detail below
 
1. Download (or update) this repository and the prerequisites.
2. Compile the code.
3. Pre-process the data.
4. Configure the model's settings.
5. Run the model.
6. Post-process the outputs.

When things don't work complain to Paul (if its a legacy bug) or Josh (if its a new bug or about this documentation).


## Prerequisites and download

1. Install relevant R packages: `install.packages(c('tidyverse', 'codas', 'assertr', 'rmarkdown'))`
2. Install the GSL and OpenMP headers. For Ubuntu this can be done with `sudo apt install libgsl-dev libomp-dev`.
3. Clone this repo


## Compilation

Normally, compile using make `make rtm_optim` in the root of the project.
If you need to specify a library location using `-L` then this can be done with `LDFLAGS` (ie: `LDFLAGS="-L<path_to_library>" make rtm optim`).
You can make `LDFLAGS` an environment variable to avoid having to do this each time.

Use the `-j<n>` option to make to specify using `n` threads.


### Further details

Three make targets are available: `rtm_optim` which compiles optimised, multi-threaded code; `rtm` which compiles unoptimised, multi-threaded code; and `rtm_debug` which compiles unoptimised, single-threaded code. For use (ie: not debugging), the `rtm_optim` target is recommended.
The standard `all` and `clean` targets are also available.



## Pre-process data

One or two data streams are currently used in the code: confirmed cases ("linelist") data and deaths data.
The current default configuration uses only deaths data.
For legacy reasons, the linelist data is called GP data within the code and the deaths data is called hospitalisation (or hosp).

1. (optional) Copy the data files to `data/raw/`. Scripts can run with data anywhere though. The data files need to be in csv format.
2. Modify `R/data/format_linelist.R` and `R/data/format_deaths.R` as required.
	* Default is that we refer to data published yesterday, change the `date.data` variable (in both files) if needed.
	* Default is that the location of the data is a command-line argument, it can also be set in the script (eg: if running from RStudio).
	* Default is ignoring 2 days of data due to reporting delay. This should be confirmed by overlaying the previous day's data and updated if required.
	* Default is running as one region covering England. `format_deaths.R` needs changing if this is not true (how?).
3. Run `Rscript R/data/format_linelist.R <path_to_data>`.
4. Run `Rscript R/data/format_deaths.R <path_to_data>`.
5. Check the output from these scripts, they warn on some obvious issues with the data. The formatted data is saved under `data/RTM_format/`

Common issues with reading the data are column names or date formats changing.
Column names can be changed using the `col.names` list.
Date formatting is inferred in the `R/data/utils.R` file, adding extra formats to the vector passed to `parse_date_time`.
Invalid rows in the death data are fixed in the `fix_dates` function.


### More detail

* Currently ``format_linelist.R`` is currently set up to format data for the United Kingdom as a whole. Output files will carry the date stamp. There will also be a trivial denominator file which sets the population from which it is possible to appear in the dataset.
* ``format_linelist.R`` also carries some code that will estimate the distribution of times from date of onset to date of lab reporting. You may wish to use this to specify the delay distribution in the code input files, discussed later. This bit of code will undoubtedly fail due to the presence of negative survival times - these are entirely possible (diagnosis while sub-clinical). Currently, I manually set these intervals to be zero.


## Configure the model

1. Check the top of the `set_up_inputs.R` file to configure most settings.
2. If you wish to change fixed parameters in the model or the priors used then edit `mod_pars.R`.
3. Run `Rscript set_up.R` to generate the config files that will actually be read by the model (in `model_runs/<run_name>`).


### Input file details

* This file will create an output directory in which the code will be run. This will be in the directory specified by ``out.dir``.
* Use ``gp.flag`` and ``hosp.flag`` to specify which data inputs are being used. Confusingly, currently I am using ``gp.flag`` to indicate if the line listing data is being used, and ``hosp.flag`` to indicate if the deaths data is being used.
* The above variables then have ``start`` and ``end`` variables indicating the range of numbered days over which these data items should be included in the likelihood. Unfortunately the range has to be continuous.
* ``Age.grps`` allows you to specify the age groups being considered, code will have to be written here when this is anything other than ``All``.
* ``ndays`` is the number of days for the run. This should include the number of days covered by the data, plus a suitable projection period.
* The contact model variables, those prefixed ``cm.`` indicate the timing of any changepoints in the contact patterns (e.g. such as due to a school holiday or a social distancing intervention), and the matrices that should be used in the periods that straddle the breakpoints. The breakpoint gives the last day that uses the current contact matrix, prior to changing. You need to specify the contact matrix and a multiplier matrix, which gives the parameterisation. See the PNAS paper for how these work - the multiplier matrix just gives a matrix of (base 0) parameter indices for the contact model parameter which will multiply the corresponding component of the contact matrix at that time.  * This file will create an output directory in which the code will be run. This will be in the directory specified by ``out.dir``.
* Use ``gp.flag`` and ``hosp.flag`` to specify which data inputs are being used. Confusingly, currently I am using ``gp.flag`` to indicate if the line listing data is being used, and ``hosp.flag`` to indicate if the deaths data is being used.
* The above variables then have ``start`` and ``end`` variables indicating the range of numbered days over which these data items should be included in the likelihood. Unfortunately the range has to be continuous.
* ``Age.grps`` allows you to specify the age groups being considered, code will have to be written here when this is anything other than ``All``.
* ``ndays`` is the number of days for the run. This should include the number of days covered by the data, plus a suitable projection period.
* The contact model variables, those prefixed ``cm.`` indicate the timing of any changepoints in the contact patterns (e.g. such as due to a school holiday or a social distancing intervention), and the matrices that should be used in the periods that straddle the breakpoints. The breakpoint gives the last day that uses the current contact matrix, prior to changing. You need to specify the contact matrix and a multiplier matrix, which gives the parameterisation. See the PNAS paper for how these work - the multiplier matrix just gives a matrix of (base 0) parameter indices for the contact model parameter which will multiply the corresponding component of the contact matrix at that time.


### Pars file details

* Specify either fixed values for certain parameters, or their initial conditions and prior parameters.
* If a new parameter is to be estimated, or conversely, is to be held fixed, this might need some editing of the ``mod_pars.Rmd`` file contained in ``./inputs/mod_pars.Rmd``.
* Refer to presentation for more details on how to specify parameters in ``mod_pars`` files.


## Run the model

1. `cd model_runs/<run_name>` (this will be fixed so that the model can run anywhere soon(TM)).
2. `../../rtm_optim` (or your version of choice, see compilation section of this document).
3. Ignore the warning about the missing `rtm_input_files.txt`.
4. Check the `adaptive_report` and `posterior_report` files to see that sensible things are happening (these will be produced periodically throughout the run).


## Post-process the outputs

1. Either `cd ../../` (ie: move to repository root) then `Rscript R/output/projections.R`, or `Rscript ../../R/output/projections.R`.
2. Outputs are PDF and RData files in `model_runs/<run_name>`.

These scripts read in values from the set up files above.
Therefore, they should not need modifying.
Blame Josh if they do....
