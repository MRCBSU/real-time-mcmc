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
<br>



[PHE gitlab]: https://gitlab.phe.gov.uk/Paul.Birrell/real-time-mcmc