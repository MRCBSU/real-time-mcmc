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
0. Bring over some chains from a previous model run, to take an initial value. Typically this will be an `.RData` file that contains an MCMC object limited to only the last iteration of the chain (to limit file size).
1. (If there is a new supply of contact matrices) Copy the .zip file containing the matrices into a new sub-folder within ``contact_mats``. The variable ``matrix.dir`` in the file ``set_up_inputs.R`` will have to be set to this name. Currently this is changed from week to week by updating the date held by the variable ``google.data.date``. A base_matrices.rds file needs to be copied into or linked to from within this directory. The most up-to-date .rds file should be stored in ``contact_mats/base_matrices``. By default, contact matrices will not overwrite existing matrices with the same file name.
2. (If there is a new supply of prevalence data) Save (or link to) the outputs from the ONS SRS in the folder ``data/raw/prevalence``. The code will pick up the most recent data file contained in this directory, however the filename has to have the format ``{date}{geography}_data.csv``, where {date} is replaced by a yyyymmdd date and {geography} takes values either ``nhs`` or ``ons`` depending on how the country is to be divided. Update the variable ``num.prev.days`` to give the length of the prevalence dataset. It is recommended to set the variable ``date.prev`` in ``config.R`` to the date on this data file.
3. If using data on deaths, go to 3a. If using data on hospital admissions, go to 3b.
   a. Set ``hosp.flag <- 1``. Set ``data.desc <- "deaths"``. Data on deaths are fetched from the modelling cell, where the death line list is usually saved as an ``.xlsx`` file. Copy (or link to) a ``.csv`` version of this file in the directory ``data/raw/deaths`` within the repo. When the data processing code runs, it will flag up some data rows that are being omitted due to dates of onset that are inconsistent with dates of deaths. Based on the output produced when doing this, I go back and clean some of the data in some of the rows suggested (usually where the date of symptom onset is more than a week beyond the date of death are these are usually due to easily-identifed typos). On line 11, I set ``args`` to be equal to ``today() - days(n)`` where ``n`` is the number of days prior to today that the deaths line list being used was generated.
   b. Set ``hosp.flag <- 0``. Set ``data.desc <- "admissions"``. Data on admissions are downloaded from Seb Funk's website at LSHTM, downloaded as an `.rds` file and then renamed to give it a date equal to the latest date in the data +1 day. This file should be saved, or linked to, from the directory ``data/raw/admissions`` within the repo. On line 11, I set ``args`` to be equal to ``today() - days(n)`` where ``n`` is the number of days prior to today that the admissions dataset is dated.
2. Serological data - These data will typically be lifted straight from the JMT shared area (or from whichever directory is specified by the `input.loc' variable defined at the top of the `format.sero.R` file. There are a number of variables whose values might need changing in `config.R`: `NHSBT.flag` typically will be set to 1 to use the NHSBT serology collection from the second and third waves. `RocheS.flag` defines which assay gets used in the data for the second and third wave data; `sero.end.date` the last possible date from which serology data would be used; `sero.date.fmt` defines the format of dates in the data file; `fix.sero.test.spec.sens` a flag indicating whether we want to use fixed estimates for serological sensitivity and specificity, or whether these should be estimated.
5. Vaccination data will be downloaded automatically from the modelling cell directory. These data take a considerable time to download as the file is large (many millions of rows). Therefore in ``config.R`` one can set ``vac.overwrite`` to ``FALSE`` if the data import routines have been previously run on the current data so that they only need to be downloaded once. If projections are required, then the variable ``future.n`` in ``config.R`` should be set to the number of vaccines anticipated be given over the coming few weeks. If ``str.date.vacc`` exists the code will expect vaccination data from this date to have been pre-prepared and formatted and to be ready available already within the repo.

## Configure the model

1. Almost all parameters that will change between runs are set in ``config.R``. No changes other than those already mentioned in the above paragraphs should be necessary if simply extending a single analysis.

### Input file details

* This file will create an output directory in which the code will be run. This will be in the directory specified by ``out.dir``.
* Use ``gp.flag`` and ``hosp.flag`` to specify which data inputs are being used. Confusingly, currently I am using ``gp.flag`` to indicate if the line listing data is being used, and ``hosp.flag`` to indicate if the deaths data is being used.
* The above variables then have ``start`` and ``end`` variables indicating the range of numbered days over which these data items should be included in the likelihood. Unfortunately the range has to be continuous.
* ``Age.grps`` allows you to specify the age groups being considered, code will have to be written here when this is anything other than ``All``.
* ``ndays`` is the number of days for the run. This should include the number of days covered by the data, plus a suitable projection period.
* The contact model variables, those prefixed ``cm.`` indicate the timing of any changepoints in the contact patterns (e.g. such as due to a school holiday or a social distancing intervention), and the matrices that should be used in the periods that straddle the breakpoints. The breakpoint gives the last day that uses the current contact matrix, prior to changing. You need to specify the contact matrix and a multiplier matrix, which gives the parameterisation. See the PNAS paper for how these work - the multiplier matrix just gives a matrix of (base 0) parameter indices for the contact model parameter which will multiply the corresponding component of the contact matrix at that time.
* This file will create an output directory in which the code will be run. This will be in the directory specified by ``out.dir``.
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

1. Change directory to `model_runs/<yyyymmdd>' where `yyyymmdd' is the data of the set of model runs. Run `R/output/chain_comparison.R' for each pair of chains. This will produce some superficial convergence checks, generating .pdf plots in the folder corresponding to the first chain. To save disk space, I'm going to presume that one chain is chosen from the two as the better/more converged.
2. For the selected chains, run the post-processing to generate the summary .html from within the `yyyymmdd' directory. This involves running the script `../../submission_scripts/submit_post_process'. In the preamble, the `array' variable needs to be changed to select only the directories that are the most converged (the directories form a base-0 list in the order in which they appear when using the `ls' command).
3. Simulate from the model. Typically this will need to be done weekly from ONS60, NHS60 and NHS28 runs. Within the `yyyymmdd' directory, edit the `../../submission_scripts/submit_simulate' so that the `array' variable looks at the right directories, and edit the variable `forecast_switch' to be equal to ``midterm''. This script will call an `R' file which also needs to be checked. Within `../../R/output/mc_midterm_forecast.R' the variable `projections.basename' should be set to "projections_midterm".
4. For the ONS60 we needs a projections report. Again, from the `yyyymmdd' directory, need to run the script `../../R/output/submit_projections_report' which should only require editing to correct the array number of the desired directory.
5. For the NHS60 we need to calculate the ICU bed occupancy plots. Download the latest NHS Daily Covid SitRep report. The data that are required are in the sheet `Org level data (inc ISP)'. Look at the cells corresponding to the NHS regions under the columns labelled `Total number of beds, as at 08:00 (occupied with confirmed COVID patients) (Total)'. Open the file `../../R/output/hosp_thresholds.R' and edit the variables
   * `output.dir' to give the location of the NHS60 run to be used
   * `beds.used' to give the contents of the column described above
   * `beds.day' to be the date of the day before the data on the SitRep document
   Source the file into R and it will produce six `.png' plots. These are usually e-mailed to Nick Gent and Duncan Kerrod on Monday mornings alongside the two .html reports generated for the `ONS60' run.
6. Navigate to the ONS60 output directory and run the submission script `../../../submission_scripts/submit_phe'. This will generate the data for the PHE prevalence report which is to go into the JMT shared area in the `Real time model outputs/data' folder.
7. Run the counterfactual simulations from the ONS60 analysis. As above for the use of the `submit_simulate' script, but in the file `mc_midterm_forecast.R' change the variable `projections.basename' to "projections_counter". Run the `submit_simulate' script from the same location as above with the correct array number.
8. Once the counterfactual simulations are complete, navigate to the directory in which the simulations have been run and source the R script `../../../R/output/vac.comparison.R'. This will (currently) generate an .RData file that is to go to the JMT shared area.
9. For SPI-M, navigate to the directory in which the model outputs exist. Open `../../../R/output/spim_format.R'. Make sure ``med.term.flag`` and ``nowcast.flag`` are both set to ``TRUE``. Update the variable ``mtp.filter.date`` by adding one week to the date (typically this should be close to ten days prior to the current date).
