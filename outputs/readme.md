```country_data.rds``` is a list that contains country-specific data used in the study. It consists of 
1) A matrix containing data for all countries included in the epidemic simulations 
2) A matrix containing country-specific data on sequencing rates as estimated from GISAID 
3) A list containing the distribution of turnaround times for each country in the epidemic simulations, as estimated from GISAID submissions
4) The mobility matrix used in the epidemic simulations, giving the daily per-capita rate of travel between each pair of countries

```time_to_detection.rds``` contains output for the single-country analyses presented in Fig. 2.

```effective_increase.rds``` contains output of the analysis that quantifies the equivalent fold increase for each reduction in turnaround time.

```GLEAM_validation.rds``` contains data for the validation of the epidemic model against GLEAM.
The folder ```GLEAM_validation``` contains raw outputs for the GLEAM simulations used to validate the metapopulation model.

```combined_outputs.rds``` contains outputs for global genomic surveillance simulations for variant detection.

```lead_times.rds``` contains genomic simulation output for the time between a variant’s first global detection and its arrival in a focal country.



