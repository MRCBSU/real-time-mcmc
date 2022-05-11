## How to merge code to drop burn in section for reports and other generated outputs

To merge this code it can be done through either merging the branch *preprocess_clean* or *SMED_Test_JK_Branch*. The rest of this explanation assumes that the reader merges using the former branch.

### Files to merge

The following files should be merged:

- All (new) files in {repo_root}/R/postprocess:
    - *create_rdatas.R*
        - Creates the Rdata which have dropped all of the burn in iterations
    - *create_report.R*
        - This is a replacement for the PostProcess.R postprocessing file in the root of the repository.
    - *traces_mcmc.R*

- Additional (updates to) files in {repo_root}/R/output:
    - *prev_report.Rmd*
    - *report_updated.Rmd*
    - *results_api.R*
    - *tidy_output.R*

### Further modifications:

The following modifications need to be made:

- Modifications to {repo_root}/R/postprocess/create_report.R
    ```r
    # Line 25
    file.loc <- dirname(thisFile())

    # Lines 29-31
    if (!file.exists("mcmc.RData")) {
        source(file.path(Rfile.loc, "tracePlots.R"))
        source(file.path(dirname(Rfile.loc), "postprocess", "create_rdatas.R"))
    }
    ```
- Modifications to {repo_root}/R/output/results_api.R
    ```r
    # Line 28
    load(file.path(out.dir, "mcmc.RData"))
    ```
-  Modifications to {repo_root}/PostProcess.R
    ```r
    # Line 19 (replaces rerun.R)
    source(file.path(dirname(thisFile()), "R", "postprocess", "create_report.R"))
    ```