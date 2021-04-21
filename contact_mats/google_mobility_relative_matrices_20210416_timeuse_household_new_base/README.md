# Scripts contained within

## For SPI-M submissions

### Creating submissions

* **spim_aggregate.R**: edit this file with the locations of results for the four nations, then it will output a combined CSV. Need to run from the 'real-time-mcmc' folder, although the results can be anywhere as long as entered in this file. For the non-English nations, the MTP are assumed to be contained in the mcmc.RData file; for England it is assumed to be in projections_midterm.RData. These assumptions should easily be changeable.
	* To use: first edit file, then `cd <path>/real-time-mcmc`, and then `Rscript <path>/spim_aggregate.R`.
* **spim_aggregate_projections.R**: as previous but use projections RData files for everywhere.
* **spim_combine.R**: take an arbitrary number of SPI-M formatted CSV files and outputs a combined CSV.
	* Usage: `Rscript <path>/spim_combine.R <in_file> <in_file> ... <out_file>`
* **spim_format.R**: output current working directory in SPI-M CSV format. Outputs nowcast and MTP assuming all contained in mcmc.RData.
	* Usage: `cd real-time-mcmc/model-runs/<path> && Rscript spim_format.R`.
* **format_midterm.R**: output MTP from current working directory in SPI-M CSV format. Outputs assuming all contained in projections_midterm.RData.
	* Usage: `cd real-time-mcmc/model-runs/<path> && Rscript format_midterm.R`.
* **spim_projections_format.R**: as previous, takes two arguments. First is the RData file, second is name of CSV file to output.
* **spim_projections_format.R**: as previous, takes two arguments. First is the RData file, second is name of CSV file to output.

### Validator issues

* **correct_DSTL_validator_errors.R**: corrects the most common validator errors in in_file and outputs a new CSV as out_file.
	* Usage: `Rscript <path>/correct_DSTL_validator_errors.R <in_file.csv> <out_file.csv>`
* **add_scenario.R**: add or change the scenario field for a whole file (NB: overwrites file)
	* Usage: `Rscript <path>/add_scenario.R <file.csv>`

## Moving files

* **move_contact_matrices.sh**: move contact matrices from Edwin from the current directory to the correct subfolder in `real-time-mcmc` with symlinks etc. Needs editing before each usage for the correct output location.
* **hpc_transfer.sh**: move aggregated data files in RTM format to the hpc.
* **rdata_to_hpc.sh** and **rdata_to_morricone.sh**: move a specific runs RData files to either the hpc or morricone for further processing (eg: to do scenarios on the DAs on the cluster).
* **send_to_morricone.sh**: use rsync to send the important files in `model_runs` to the projects folder on morricone.

