# real-time-mcmc

Paul's current process for updating model results

## Process overview

Quick summary of basic steps
 
1. Download (or update) this repository and the prerequisites.
2. Compile the code.
3. Pre-process the data.
4. Configure the model's settings.
5. Run the model.
6. Post-process the outputs.

## Pre-process data

Five data streams are currently used in the code, with the option to include a fifth (confirmed cases).
- Counts of deaths
- Serological testing data
- Prevalence estimates
- Counts of vaccinations
- Contact matrices provided by Edwin van Leeuwen's analysis of Google mobility data, DfE attendance data and UKTUS data.
For legacy reasons, the linelist data is called GP data within the code and the deaths data is called hospitalisation (or hosp).

Re-running the model isn't quite as simple as just running the model code (which will seek out the latest data). Model set-up and data processing is typically done in a secure environment outputs from which are then transferred to the HPC.
1. (If there is a new supply of contact matrices) Copy the .zip file containing the matrices into a new sub-folder within ``contact_mats``. The variable ``matrix.dir`` in the file ``set_up_inputs.R`` will have to be set to this name. Currently this is changed from week to week by updating the date held by the variable ``google.data.date``. A base_matrices.rds file needs to be copied into or linked to from within this directory. The most up-to-date .rds file should be stored in ``contact_mats/base_matrices``.
   2. (If there is a new supply of prevalence data) Save (or link to) the outputs from the ONS SRS in the folder ``data/raw/prevalence``. The code will pick up the most recent data file contained in this directory, however the filename has to have the format ``{date}_{days}days_{geography}_data.csv``, where {date} is replaced by a yyyymmdd date, {days} gives the length in days of the ONS output and {geography} takes values either ``nhs`` or ``ons`` depending on how the country is to be divided. Update the variable ``num.prev.days`` to give the length of the prevalence dataset. It is recommended to set the variable ``date.prev`` in ``config.R`` to the data on this data file.
3. Data on deaths is fetched from the modelling cell, where the death line list is usually saved as an ``.xlsx`` file. Copy (or link to) a ``.csv`` version of this file in the directory ``data/raw/deaths`` within the repo. When the data processing code runs, it will flag up some data rows that are being omitted due to dates of onset that are inconsistent with dates of deaths. Based on the output produced when doing this, I go back and clean some of the data in some of the rows suggested (usually where the date of symptom onset is more than a week beyond the date of death are these are usually due to easily-identifed typos). On line 11, I set ``args`` to be equal to ``today() - days(n)`` where ``n`` is the number of days prior to today that the deaths line list being used was generated.
2. Serological data remain untouched for now. Nothing to change here.
5. Vaccination data will be downloaded automatically from the modelling cell directory. These data take a considerable time to download as the file is large (many millions of rows). Therefore in ``config.R`` one can set ``vac.overwrite`` to ``FALSE`` if the data import routines have been previously run on the current data so that they only need to be downloaded once. If projections are required, then the variable ``future.n`` in ``config.R`` should be set to the number of vaccines anticipated be given over the coming few weeks.

## Configure the model

1. Almost all parameters that will change between runs are set in ``config.R``. No changes other than those already mentioned in the above paragraphs should be necessary if simply extending a single analysis.

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
