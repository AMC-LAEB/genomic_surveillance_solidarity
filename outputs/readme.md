**simulation_onsets.xlsx** contains each country's day of first variant infection, for each of the 10,000 simulations in the main text, for each value of variant R<sub>e</sub>.

**country_data.rds** is a list that contains country-specific data used in the study. It consists of 
1) A matrix containing data for all countries included in the epidemic simulations (columns: country, population size, continent, centroid longitude, centroid latitude, sequencing rate)
2) A matrix containing country-specific data on sequencing rates as estimated from GISAID (columns: country, population size, continent, SARS-CoV-2 sequence count, SARS-CoV-2 median turnaround time, SARS-COV-2 sequencing rate, per capita GDP, World Bank income classification, influenza virus sequence count, influenza virus sequencing rate)
3) A list containing the distribution of turnaround times for each country in the epidemic simulations, as estimated from GISAID submissions
4) The mobility matrix used in the epidemic simulations, giving the daily per-capita rate of travel between each pair of countries

**detection_outputs_threshold_F.rds** contains outputs for each the global genomic surveillance simulations for variant detection.

**detection_outputs_threshold_T.rds** contains outputs for each the global genomic surveillance simulations for the time until the variant is found to exceed a variant proportion threshold.

**effective_increase.rds** contains output of the analysis that quantifies the equivalent fold increase for each reduction in turnaround time.

**GLEAM_validation.rds** contains data for the validation of the epidemic model against GLEAM.

**lead_times.rds** contains genomic simulation output for the time between a variant’s first global detection and its arrival in a focal country.

**time_to_detection.rds** contains output for the single-country analyses presented in Fig. 2.

The folder **GLEAM_validation** contains the outputs for the GLEAM simulations used to validate the metapopulation model.
