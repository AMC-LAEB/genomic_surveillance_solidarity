## Running the analyses

Analyses require the R packages ```data.table```, ```tidyr```, ```stringr```, ```Rcpp```, ```doParallel``` and ```npreg```, as well as ```GNU Parallel```. Note: please ensure your R session is running in an English locale, e.g. by running ```Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")```.

1. Run ```getCountryData.R``` to analyze the global SARS-CoV-2 sequence metadata. This requires a file containing metadata for global sequence data, which can be extracted from GISAID. The output file ```country_data.rds``` can be found in the ```../outputs``` folder and can be placed in the current folder to reproduce the further analyses.

2. Run ```singleCountrySimulations.R``` to generate the outputs for the single-country simulations presented in Figure 2.

3. Run ```simulateGlobalEpidemics.R``` to generate the global metapopulation epidemic simulations. Note: by default, these simulations are set to be run on 50 cores.

4. Run ```surveillanceSimulations.bash``` to simulate the global genomic surveillance process. This runs ```runSurveillanceSimulations.R``` in parallel for the different surveillance parameters, with the core of the simulation performed in ```surveillanceSimulations.cpp```, followed by ```getLeadTimes.R``` to compute the time between global detection and local arrival. The raw epidemic simulation output used as input can be downloaded from https://zenodo.org/records/10051237 and should be put in a folder ```simulations```. Note: by default, these simulations are set to be run on 80 cores.

Specific analyses corresponding to individual figures can be found in the ```../figures``` folder. 


