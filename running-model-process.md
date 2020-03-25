# Steps to go through when running the model during the CoVID-19 outbreak.

To get the RTM package, make a branch of the software package currently available from [PHE gitlab].

The below steps do not have to be adhered to, rather they are my current workflow for iterating the analyses as new data become available each day.

* Receive new datasets via e-mail. Save the files to the relevant sub-directory of the `data streams' directory on the Modelling cell shared workspace [//phe.gov.uk/porton/SharedData/ERD/Modelling\ cell\ nCoV\ 2019], if not already saved to there.



* From the original data files, I save local .csv copies of the sheets containing the data needed.
  * To minimise editing of input routines, the local .csv file for the death data should be named '20200324 COVID19 Deaths.csv', where the 8-digit date code changes each day.
  * To minimise editing of input routines, the local .csv file for the line list data should be named '20200324 - Anonymised Line List.csv', where the 8-digit date code changes each day.
<br>

* In the repository, go to R/data to find the data formatting R routines. There are three types of variable that need to be checked before running any of these files and these are printed at the head of the .R files.
  * The top are inputs that should change on a daily basis (i.e. the date).
  * The second will be variables that are dependent on the user (i.e. where the data has been saved locally).
  * The third are variables that, unfortunately, commonly change in the input data files, should as field names and date formats. These need to be checked against the input file.
<br>


* Bearing in mind the above points, run format_deaths.R and format_linelist.R. These will save data outputs in the repository in './data/deaths/' and in './data/Linelist/'.
  * Common bugs that may arise usually centre on the way in which the dates are being read in from the input (.csv) data file. This can cause functions such as ``lubridate::as_date()`` to return ``NA`` values. In this event, it might be necessary to adapt the code to convert them to character string and then back again to dates.
  * Currently ``format_linelist.R`` is currently set up to format data for the United Kingdom as a whole. Output files will carry the date stamp. There will also be a trivial denominator file which sets the population from which it is possible to appear in the dataset.
  * ``format_linelist.R`` also carries some code that will estimate the distribution of times from date of onset to date of lab reporting. You may wish to use this to specify the delay distribution in the code input files, discussed later. This bit of code will undoubtedly fail due to the presence of negative survival times - these are entirely possible (diagnosis while sub-clinical). Currently, I manually set these intervals to be zero.
  * In ``format_deaths.R`` be careful to check the matching of the region and date of death fields.
<br>

* Return to the root of the repository.

* Open up the R files beginning `set_up*`.

* In ``set_up.R``:
  * Change the date of the runs to the date being used for the data.
  * List the regions in which you're running the model. These much match to a geographical unit specified in ``set_up_inputs.R``. These will be typically the same as in the ONS population data (contained within the repository). Country-level units are typically written in capitals. Use underscores for spaces.
  * Use these in the variable ``out.dir`` to give a name to the run.
  * If you want to run the code from within R, uncomment the system command, the penultimate command in the file.

* In ``set_up_inputs.R``:
  * This file will create an output directory in which the code will be run. This will be in the directory specified by ``out.dir``.
  * Use ``gp.flag`` and ``hosp.flag`` to specify which data inputs are being used. Confusingly, currently I am using ``gp.flag`` to indicate if the line listing data is being used, and ``hosp.flag`` to indicate if the deaths data is being used.
  * The above variables then have ``start`` and ``end`` variables indicating the range of numbered days over which these data items should be included in the likelihood. Unfortunately the range has to be continuous.
  * ``Age.grps`` allows you to specify the age groups being considered, code will have to be written here when this is anything other than ``All``.
  * ``ndays`` is the number of days for the run. This should include the number of days covered by the data, plus a suitable projection period.
  * The contact model variables, those prefixed ``cm.`` indicate the timing of any changepoints in the contact patterns (e.g. such as due to a school holiday or a social distancing intervention), and the matrices that should be used in the periods that straddle the breakpoints. The breakpoint gives the last day that uses the current contact matrix, prior to changing. You need to specify the contact matrix and a multiplier matrix, which gives the parameterisation. See the PNAS paper for how these work - the multiplier matrix just gives a matrix of (base 0) parameter indices for the contact model parameter which will multiply the corresponding component of the contact matrix at that time.
<br>

* In ``set_up_pars.R``:
  * Specify either fixed values for certain parameters, or their initial conditions and prior parameters.
  * If a new parameter is to be estimated, or conversely, is to be held fixed, this might need some editing of the ``mod_pars.Rmd`` file contained in ``./inputs/mod_pars.Rmd``.
  * Refer to presentation for more details on how to specify parameters in ``mod_pars`` files.
<br>

* Run the model either from the command line, via batch submission script, or even from within R.

* Complain about my coding when you immediately hit a bug.
  

[PHE gitlab]: https://gitlab.phe.gov.uk/Paul.Birrell/real-time-mcmc