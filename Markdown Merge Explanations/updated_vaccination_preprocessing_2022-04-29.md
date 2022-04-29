## Changes to vaccination preprocessing to keep in line with the increase in the number of vaccinations being delivered.

Presently it is unknown what the proportion of vaccinations delivered to the 5-14 age group will converge to. Due to a new policy to increase vaccinations in those under 12 years old and an uptick seen in the number of vaccinations delivered to this age group seen in the vaccinations linelist, it has been decided to increase the number of children expected to receive a dose by the model.

The changes to the code which reflect this are shown below:

```diff
diff --git a/Process_Vaccinations.R b/Process_Vaccinations.R
index 08a24252c3d02e40650362841c920db7a06e4750..6774ab1cb2f2ebe9fd8c9b55c2b96829c9e8370f 100644
--- a/Process_Vaccinations.R
+++ b/Process_Vaccinations.R
@@ -67,13 +67,13 @@ if (id == 0) {
     vacc.over50s <- 0.87
     vacc.25.44 <- 0.72
     vacc.15.24 <- 0.66
-    vacc.5.14 <- 0.16
+    vacc.5.14 <- 0.30
```

This value may change over the next few weeks as we observe the trend of those receiving a vaccine. 

The values were chosen to reflect the increase in the proportion vaccinated in this group (according to the line list and ONS population estimates) from 19.97% on the 5th April 2022 to 23.47% on the 27th of April 2022.