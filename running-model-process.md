# Steps to go through when running the model during the CoVID-19 outbreak.

To get the RTM package, make a branch of the software package currently available from [PHE gitlab].

The below steps do not have to be adhered to, rather they are my current workflow for iterating the analyses as new data become available each day.

* Receive new datasets via e-mail. Save the files to the relevant sub-directory of the `data streams' directory on the Modelling cell shared workspace [//phe.gov.uk/porton/SharedData/ERD/Modelling\ cell\ nCoV\ 2019], if not already saved to there.
  * To minimise editing of input routines, the local .csv file for the death data should be named '20200324 COVID19 Deaths.csv', where the 8-digit date code changes each day.
  * To minimise editing of input routines, the local .csv file for the line list data should be named '20200324 - Anonymised Line List.csv', where the 8-digit date code changes each day.

* From the original data files, I save local .csv copies of the sheets containing the data needed.

[PHE gitlab]: https://gitlab.phe.gov.uk/Paul.Birrell/real-time-mcmc