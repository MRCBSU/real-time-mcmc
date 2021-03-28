# real-time-mcmc

Josh's current process for updating model results.
For lots of detail on the different parts, see paul_workflow.md; here, I just list the standard steps I take on a Friday.

## Setup repo (on lutyens)

1. First, I make sure that any changes are fully updated in the COVID branch from lutyens, hpc, or anywhere else.
2. Checkout the latest COVID branch, and then take a run-specific branch.
	1. `git checkout COVID`
	1. `git pull`
	1. `git checkout -b <run_name>` (eg: josh_20210326)

## Update the data sources (on lutyens)

### Deaths and serology

Deaths should be available on lutyens at `/data/covid-19/data-raw/deaths/2021-03-26.csv`.
Peter transfers these each day at about 7pm.
This file is picked up automatically by the relevant R scripts so you shouldn't need to do anything.

The serology is only rarely updated, and again should be picked up automatically from the correct location on lutyens.

### Prevalence

1. Download to your local PC the csv file extracted from the ONS SRS (normally emailed by Paul or Josh on Thursday or Friday).
2. Copy (via scp) this to lutyens, the file path should be `<root_dir>/data/raw/prevalence/ONS_regions/INLA_<date>.csv` where `root_dir` is the location of your real-time-mcmc folder and `date` is the date of the last data in the file (normally the Monday). If the data is from the old MCMC model, omit "INLA_" in the filename. Eg: `~/real-time-mcmc/data/raw/prevalence/ONS_regions/INLA_2021-03-22.csv`.
3. Update the `date.prev` variable in `config.R` to the date the file is labelled with.

### Vacciantions

1. Find the latest file labelled "immunisations" in the subfolders under `/data/covid-19/data-raw/dstl/` (I've not quite worked out the reporting pattern of this yet, this week it was on the Wednesday).
2. Unzip the file, keeping it in the same directory.
3. Update the `str.date.vacc` variable in `config.R` to match the date of the data.
4. Enter the vaccination projection numbers as `future.n` in `config.R`.

### Contact matrices

1. Extract the `contact_mat.tar.gz` file from the zip named `timeuse_household.zip` attached to Edwin's email (sent each Friday).
2. Copy the file to an empty directory on lutyens (I use `~/Downloads/contact_mats/`.
4. Update the date on the `OUTPUT` line in the `move_contact_matrices.sh` to the date the matrices were produced. The `MAIN_DIR` variable may also need updating to point to where your `real-time-mcmc` directory is located.
4. In a lutyens terminal, `cd` to the directory containing just `contact_mat.tar.gz`.
5. Run the `move_contact_matrices.sh` script.
6. Update the `google.data.date` variable in `config.R` to the date the matrices are from.

### Prepare formatted data and transfer to HPC

1. Commit all the above changes to `config.R` and the new contact matrices (ie: run `git add .` and `git commit -m"Update data and contact matrices`, or some similar message).
2. Ensure that `prev.flag` is set to `1` in `config.R` and `format.inputs` is set to `TRUE` in `run.R`.
3. If any other changes are necessary in `config.R` which are relevant to how the data is used, do so.
4. Run `Rscript run.R <date>` with the appropriate date to use for the deaths data, eg: `20210326`.
5. If necessary, repeat the last two steps for each different "type" of data used (eg: different timing of prevalence data, different types of deaths). You can ignore any errors regarding being unable to read a `tmp.RData` file that does not exist.
7. Run `hpc_transfer.sh` (this file may need updating if your HPC directories are setup differently to mine) to copy the formatted data to the HPC.
8. Check if you've changed any relevant files using `git status` and `git diff`, commit these changes if necessary.
9. Push your changes to the git repo (ie: `git push`).

## Running the model (on HPC)

1. `cd` to your `real-time-mcmc` directory on the HPC.
2. Update any changes made to the GitLab repo using `git fetch`.
3. Checkout the newly branch created branch using `git checkout <run_name>`. Eg: `git checkout josh_20210326`.
4. Ensure `format.inputs` is set to `FALSE` in `run.R`.
5. Load the needed modules on the HPC using `module load R/3.6 pandoc`.
6. (Optional) update `previous.run.to.use` in `config.R` to a sub-directory of `model_runs` with a recent run. A sample from the posterior of the specified run is used as the starting point for the runs you are about to start.
7. For each model configuration you wish to run:
	1. Setup the model configuration in `config.R` and any other changes required (the latter is unusual).
	2. Run `Rscript run.R <date>` (eg: `Rscript run.R 20210326`).
8. Change to the correct subdirectory of `model_runs` (eg: `cd model_run/20210326`).
9. Run `ls` and check that all the model configurations you wish to run are present.
10. Update the `array` directive in `submit_multiple_runs`. These specify which configurations you want to run (alphabetically based on your working directory when you submit the job). Note it is 0-indexed. For example, if there are four model configurations and you wish to run all of them then specify `array=0-3`.
11. Submit the jobs to the HPC using `sbatch ../../submit_multiple_runs`.
12. Check your jobs are running using `squeue -u<username>` (eg: `squeue -ujbb_50`).
13. Any issues normally occur in the first minute or so of running, so it's worth checking that every job survives at least this long (you won't get any emails alerting you about errors until every job completes or fails).

## Post-processing (on HPC)

Often the first attempt at post-processing (creating the report-updated.html with the results) will fail.
I am unsure why this occurs, but re-running it normally fixes this.
I use the submission script `submit_post_process_all`.
Similarly, `sumbit_simulate_all` can be used to create the medium-term projections (MTPs).
There are various scripts named `mc_*_forecast.R` which can be adapted to simulate specific scenarios (eg: changing matrices into the future).

## Debugging

* To run `run.R` in an interactive session, `rm(list=ls()); date.data <- "20210326"; source('run.R')` is useful (update the date as required).
* Sometimes, the column names in the input files change (eg: capitalisation), the code is setup such that the `possible.col.names` variable within each `R/data/` file can be given a vector of column names that mean the same thing.

