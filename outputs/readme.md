**simulation_onsets.xlsx** contains each country's day of first variant infection, for each of the 10,000 simulations in the main text, for each value of variant R<sub>e</sub>.

**detection_output.xlsx** contains all global genomic surveillance simulation outputs from the main text (columns: simulation index, Variant R<sub>e</sub>, country of origin, country of detection, day of detection, global infections by day of detection, global surveillance strategy). 

**country_data.rds** is a list that contains country-specific data used in the study. It consists of 
1) A matrix containing data for all countries included in the epidemic simulations (columns: country, population size, continent, centroid longitude, centroid latitude, sequencing rate)
2) A matrix containing country-specific data on sequencing rates as estimated from GISAID (columns: country, population size, continent, SARS-CoV-2 sequence count, SARS-CoV-2 median turnaround time, SARS-COV-2 sequencing rate, per capita GDP, World Bank income classification, influenza virus sequence count, influenza virus sequencing rate)
3) A list containing the distribution of turnaround times for each country in the epidemic simulations, as estimated from GISAID submissions
4) The mobility matrix used in the epidemic simulations, giving the daily per-capita rate of travel between each pair of countries

**detection_outputs.rds** contains outputs for each the global genomic surveillance simulations (columns: detection day, global infections by day of detection, country of detection, country of origin, variant R<sub>e</sub>, global surveillance strategy (denoted by min. and max. sequencing rate), initial wildtype prevalence, wildtype transmission rate, turnaround time φ, global mobility rate).

**effective_increase.rds** contains output of the analysis that quantifies the equivalent fold increase for each reduction in turnaround time, for the single-country analyses presented in Fig. 2 (columns: variant R<sub>e</sub>, wildtype R<sub>e</sub>, equivalent fold increase, reduction in turnaround time, initial wildtype prevalence).

**GLEAM_validation.rds** contains data for the validation of the epidemic model against GLEAM (columns: index country, focal country, arrival time in the focal country in GLEAM, and arrival time in the focal country in the metapopulation model).

**lead_times.rds** contains genomic simulation output for the time between a variant’s first global detection and its arrival in a focal country, as shown in Fig. 4d (columns: country, country's mean time between first global detection and first local case, continent, the country’s sequencing rate, the global genomic surveillance strategy, the variant R<sub>e</sub>).

**time_to_detection.rds** contains output for the single-country analyses presented in Fig. 2 (columns: variant R<sub>e</sub>, wildtype R<sub>e</sub>, sequencing rate (S/M/wk), initial wildtype prevalence, day of detection, reduction in day of detection per S/M/wk added, infections by the day of detection, reduction in infections by the day of detection per S/M/wk added).


The folder **GLEAM_validation** contains the outputs for the GLEAM simulations used to validate the metapopulation model.
