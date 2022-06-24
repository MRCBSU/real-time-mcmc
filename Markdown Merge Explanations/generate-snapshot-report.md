## Steps to prepare snapshot report

This means that you will need to checkout the following files from this branch (SMED_Test_JK_Branch) for use:

- {repo}/R/outputs/SnapshotReport.Rmd

- {repo}/R/outputs/SnapshotReportMaps.Rmd
    
- {repo}/R/outputs/SnapshotRender.R
    
- {repo}/submission_scripts/submit_snapshot_report

Then the following steps will need to be undertaken:

- modify the submission script to use the module system and sbatch parameters for the cambridge hpc (currently it is setup for SMED)
- Get hold of maps.zip and unzip to the root of the repository (this contains the map files for the reports)
    - Contained in email from Joel.kandiah@phe.gov.uk
    - This should leave a new directory "Maps" containing the two rds files
    - Note: these files will not be used if the package sf is not installed
- If not already installed install the package "here" and attempt to install "sf" 
        - If "sf" does not install note that map_flag must be set to false in the SnapshotReportRender.R file

In order to generate the report the following steps must be taken:

- In R/output/mc_snapshot_forecast.R change the variable `snap.date' equal to the date for which you would like to report the snapshot (must be within the period spanned by the MCMC run).
- Generate the snapshot rdata file using submit simulate, setting variable `forecast_switch="snapshot"'.
- Modify the parameter str.date in SnapshotReportRender.R to match the chosen snapshot date (and to match the filename of the RData file generated in the previous step)
    - This should also match the format shown on the snapshot RData file (typically %Y-%m-%d)
- Confirm that map_flag in SnapshotReportRender.R is set to False if sf failed to install and is set to True otherwise
- Finally call the script from the directory {repo}/model_runs/{date}/
    - Note that multiple reports can be generated at the same time through the array parameter (specified similarly to other scripts in the repository)