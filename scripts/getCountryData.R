## Returns four outputs:
### 1: A matrix consisting of data for all countries in the GTM dataset
### 2: A matrix consisting of data for all countries in the GISAID metadata dataset
### 3: A list containing the distribution of turnaround times by day for all countries in the GTM dataset 
### 4: A matrix M where M_ij contains the daily proportion of individuals in country j moving to country i 
Sys.setlocale("LC_COLLATE","en_US.UTF-8")

getCountryData <- function(){
  
  require(data.table)
  require(tidyr)
  require(stringr)
  
  ### Generate list of alternative country names
  
  country_name_list = list()
  country_name_list[[1]] = c("Cape Verde", "Cabo Verde")
  country_name_list[[2]] = c("Côte d'Ivoire", "Cote d'Ivoire","Ivory Coast","Côte d’Ivoire")
  country_name_list[[3]] = c("Czech Republic","Czechia")
  country_name_list[[4]] = c("Congo, Dem. Rep.", "Democratic Republic of the Congo","DR Congo","Congo (Democratic Republic of the)","Congo, the Democatic Republic of")
  country_name_list[[5]] = c("Egypt, Arab Rep.", "Egypt")
  country_name_list[[6]] = c("Swaziland", "Eswatini")
  country_name_list[[7]] = c("Lao PDR", "Laos","Lao People's Democratic Republic","Lao, People's Democratic Republic")
  country_name_list[[8]] = c("Macedonia", "North Macedonia","Macedonia (the former Yugoslav Republic of)","Macedonia, the former Yogoslav Republic of")
  country_name_list[[9]] = c("Congo, Rep.", "Republic of the Congo","Congo")
  country_name_list[[10]] = c("Russian Federation", "Russia")
  country_name_list[[11]] = c("St. Maarten", "Sint Maarten")
  country_name_list[[12]] = c("Slovak Republic", "Slobakia")
  country_name_list[[13]] = c("Korea, Rep.", "South Korea","Korea (Republic of)","Korea, Republic of")
  country_name_list[[14]] = c("Virgin Islands (U.S.)", "U.S. Virgin Islands","United States Virgin Islands")
  country_name_list[[15]] = c("St. Lucia", "Saint Lucia")
  country_name_list[[16]] = c("St. Pierre and Miquelon", "Saint Pierre and Miquelon")
  country_name_list[[17]] = c("St. Vincent and Grenadines", "Saint Vincent and Grenadines", "Saint Vincent and the Grenadines","St. Vincent and the Grenadines")
  country_name_list[[18]] = c("St. Kitts and Nevis", "Saint Kitts and Nevis","Saint Kitts and Nevis, Federation of")
  country_name_list[[19]] = c("Slovak Republic", "Slovakia")
  country_name_list[[20]] = c("Micronesia, Fed. States of", "Micronesia","Micronesia (Federated States of)","Micronesia, Fed. Sts.")
  country_name_list[[21]] = c("Korea, Dem. Rep.", "North Korea","Korea (Democratic People's Republic of)","Korea, Dem. People's Rep.")
  country_name_list[[22]] = c("Wallis and Futuna Islands", "Wallis and Futuna")
  country_name_list[[23]] = c("Curaçao", "Curacao")
  country_name_list[[24]] = c("Faeroe Islands", "Faroe Islands")
  country_name_list[[25]] = c("São Tomé and Principe", "Sao Tome and Principe","São Tomé and Príncipe")
  country_name_list[[26]] = c("St-Barthélemy", "Saint Barthelemy")
  country_name_list[[27]] = c("Trinidad and Tobago", "Trinidad & Tobago")
  country_name_list[[28]] = c("United States of America", "United States", "USA")
  country_name_list[[29]] = c("Kyrgyzstan", "Kyrgyz Republic")
  country_name_list[[30]] = c("Iran", "Iran, Islamic Rep.","Iran (Islamic Republic of)","Iran, Islamic Republic of")
  country_name_list[[31]] = c("Gambia", "Gambia, The")
  country_name_list[[32]] = c("Turkey", "Turkiye")
  country_name_list[[33]] = c("Venezuela", "Venezuela, RB", "Venezuela (Bolivarian Republic of)","Venezuela, Bolivarian Republic of")
  country_name_list[[34]] = c("Syria", "Syrian Arab Republic")
  country_name_list[[35]] = c("Yemen", "Yemen, Rep.")
  country_name_list[[36]] = c("Bolivia", "Bolivia (Plurinational State of)","Bolivia, Plurinationial State of")
  country_name_list[[37]] = c("Brunei", "Brunei Darussalam")
  country_name_list[[38]] = c("United Kingdom", "United Kingdom of Great Britain and Northern Ireland")
  country_name_list[[39]] = c("British Virgin Islands", "Virgin Islands (British)")
  country_name_list[[40]] = c("Moldova", "Moldova (Republic of)","Moldova, Republic of")
  country_name_list[[41]] = c("Tanzania", "Tanzania, United Republic of")
  country_name_list[[42]] = c("Vietnam", "Viet Nam")
  country_name_list[[43]] = c("Bahamas","The Bahamas","Bahamas, The")
  country_name_list[[44]] = c("Palestine","Palestine, State of","Palestinian Territory")
  country_name_list[[45]] = c("Macao","Macau","Macao SAR, China")
  country_name_list[[46]] = c("Hong Kong","Hong Kong SAR (China)","Hong Kong SAR, China","Hong Kong (SAR)")
  country_name_list[[47]] = c("Saint Martin","St. Martin (French part)")
  country_name_list[[48]] = c("St. Maarten","Sint Maarten (Dutch part)")
  country_name_list[[49]] = c("Taiwan","Taiwan, China")
  
  
  # Read mobility data
  # https://data.jrc.ec.europa.eu/dataset/2ca0fa06-856a-49eb-9a14-797112cc6302
  
  gtm = read.csv("gtm.csv")
  
  # Rename countries
  for (i in 1:length(country_name_list)){ 
    whch_s = which(gtm$source_name %in% country_name_list[[i]])
    whch_t = which(gtm$target_name %in% country_name_list[[i]])
    if (length(whch_s)>0){
      gtm[whch_s,]$source_name = country_name_list[[i]][1]
      gtm[whch_t,]$target_name = country_name_list[[i]][1]
    }
  }
  
  ### Get country data
  
  mobility_countries = data.frame(Country=unique(gtm$source_name))
  
  # Read population size data
  # https://www.kaggle.com/datasets/iamsouravbanerjee/world-population-dataset
  
  pops = read.csv("world_population.csv")

  ctrnames = pops[,3]
  
  # Add population size data for countries with the same name in the population size dataset and the GTM dataset
  mobility_countries$Population = pops[,6][match(mobility_countries[,1],ctrnames)] 
  mobility_countries$Continent = pops[,5][match(mobility_countries[,1],ctrnames)]
  
  # For other countries, match the names in the two datasets
  for (i in 1:length(country_name_list)){
    country_pop = pops[,6][match(country_name_list[[i]],ctrnames)]
    country_pop = country_pop[!(is.na(country_pop))]
    country_continent = pops[,5][match(country_name_list[[i]],ctrnames)]
    country_continent = country_continent[!(is.na(country_continent))]
    if (length(country_pop)>0){
      mobility_countries$Population[match(country_name_list[[i]][1],mobility_countries[,1])] = country_pop
      mobility_countries$Continent[match(country_name_list[[i]][1],mobility_countries[,1])] = country_continent
    }
  }

  # This file, extracted from GISAID, contains the submission date, collection data and sampling locations of all viruses in the GISAID EpiCoV database
  mtdt_cov = fread("dates_and_locations.tsv",sep='\t')
  
  # Extract country of sampling
  mtdt_cov$Location = paste0(mtdt_cov$Location, " /") 
  mtdt_cov$Country = str_match(mtdt_cov$Location, "/\\s*(.*?)\\s*/")[,2]
  
  # Extract turnaround time
  mtdt_cov$`Collection date` = as.Date(mtdt_cov$`Collection date`,'%Y-%m-%d')
  mtdt_cov$`Submission date` = as.Date(mtdt_cov$`Submission date`,'%Y-%m-%d')
  mtdt_cov$TAT = difftime(mtdt_cov$`Submission date`,mtdt_cov$`Collection date`,units='days')
  
  # Keep only the year 2022
  sequencing_rate_data = as.data.frame(table(mtdt_cov[mtdt_cov$`Collection date`>="2022-01-01" & mtdt_cov$`Collection date`<"2023-01-01" & mtdt_cov$`Submission date`< "2023-07-01",]$Country))
  colnames(sequencing_rate_data) = c("Country","Count")
  sequencing_rate_data[,1] = as.character(sequencing_rate_data[,1])
  
  # Rename countries to ensure consistency
  old_names = sequencing_rate_data[,1]
  for (i in 1:length(country_name_list)){
    whch_s = which(sequencing_rate_data$Country %in% country_name_list[[i]])
    if (length(whch_s)>0){
      sequencing_rate_data[whch_s,]$Country = country_name_list[[i]][1]
    }
  }
  
  # Get turnaround time information
  
  # Initialize output
  tats_by_day = list()
  for (i in 1:length(mobility_countries[,1])){
    tats_by_day[[i]] = rep(0,500)
  }
  sequencing_rate_data$Median_TAT = NA
  
  # Get correpondence between indices in GTM dataset and GISAID dataset
  idxes = match(sequencing_rate_data[,1],mobility_countries[,1])
  
  for (i in 1:nrow(sequencing_rate_data)){
    ctry = old_names[i] # Get country name
    if (length(mtdt_cov[mtdt_cov$Country==ctry,]$TAT)>0){
      tats = mtdt_cov[mtdt_cov$Country==ctry,]$TAT # Get all turnaround time values for this country
      median_tat = median(tats,na.rm=T) # Get median value
      tat_by_day = tabulate(as.numeric(tats)+1,nbins=500)/length(tats[!(is.na(tats))]) #For each day get proportion of sequences with turnaround time equal to that day
      if (!(is.na(idxes[i]))){
        tats_by_day[[idxes[i]]] = tat_by_day # Save distribution
      } 
      sequencing_rate_data[i,]$Median_TAT = median_tat # Keep median
    }
  }

  ### Add population size data
  
  #Combine countries from GTM dataset with countries from GISAID dataset, to include
  #countries with no sequences in GISID
  
  sequencing_rate_data = merge(mobility_countries[,1:3],sequencing_rate_data,all=T)
  
  ctrnames = pops[,3] # Get country names
  
  # Add population size data for countries with the same name in the population size dataset and the GISAID dataset
  sequencing_rate_data$Population = pops[,6][match(sequencing_rate_data[,1],ctrnames)]
  sequencing_rate_data$Continent = pops[,5][match(sequencing_rate_data[,1],ctrnames)]
  
  # For other countries, match the names in the two datasets
  for (i in 1:length(country_name_list)){
    country_pop = pops[,6][match(country_name_list[[i]],ctrnames)]
    country_pop = country_pop[!(is.na(country_pop))]
    country_continent = pops[,5][match(country_name_list[[i]],ctrnames)]
    country_continent = country_continent[!(is.na(country_continent))]
    if (length(country_pop)>0){
      sequencing_rate_data$Population[match(country_name_list[[i]][1],sequencing_rate_data[,1])] = country_pop
      sequencing_rate_data$Continent[match(country_name_list[[i]][1],sequencing_rate_data[,1])] = country_continent
    }
  }
  
  ### Compute sequencing rate in S/M/wk
  sequencing_rate_data$Seqrate = sequencing_rate_data$Count/(sequencing_rate_data$Population/1e6)/52
  
  ### Get income classification
  #https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups
  
  incomes = read.csv("income_levels.csv")
  ctrnames = incomes[,1]
  sequencing_rate_data$Income =incomes[,4][match(sequencing_rate_data[,1],ctrnames)] #Add income classification data for countries with the sane name in the population size dataset and the GISAID dataset
  
  # For other countries, match the names in the two datasets
  for (i in 1:length(country_name_list)){
    country_income = incomes[,4][match(country_name_list[[i]],ctrnames)]
    country_income = country_income[!(is.na(country_income))]
    if (length(country_income)>0){
      sequencing_rate_data$Income[match(country_name_list[[i]][1],sequencing_rate_data[,1])] = country_income
    }
  }


  sequencing_rate_data[is.na(sequencing_rate_data$Count),c("Count","Seqrate")] = 0

  # Get each country in the mobility dataset's sequencing rate 
  mobility_countries$Seqrate = sequencing_rate_data$Seqrate[match(mobility_countries$Country,sequencing_rate_data$Country)]
  
  # Remove subnational countries in the mobility dataset for which the GISAID metadata is unrepresentative
  sequencing_rate_data = sequencing_rate_data[!(sequencing_rate_data$Country%in%c("Gibraltar","Macao","Niue")),]

  # Reorder mobility data
  ord = order(mobility_countries[,1])
  mobility_countries = mobility_countries[ord,]
  tats_by_day = tats_by_day[ord]
  
  ### Generate mobility matrix for use in global epidemic simulations
  
  mobility_data = gtm[gtm$year==2016,] # Most recent year
  mobility_data = mobility_data[,c(1,2,6)]
  mobility_data_wide = spread(mobility_data,key="target_name",value="estimated_trips")
  mobility_data_wide = mobility_data_wide[,2:ncol(mobility_data_wide)]
  rownames(mobility_data_wide) = colnames(mobility_data_wide) 
  mobility_data_wide =(mobility_data_wide+t(mobility_data_wide))/2 # Symmetrize
  for (i in 1:nrow(mobility_data_wide)){
    ### For each country i, compute the daily probability that an individual from country i moves to country j
    mobility_data_wide[,i] = (mobility_data_wide[,i]/mobility_countries$Population[i])/365
  }
  
  # Add country names to TAT distributions
  names(tats_by_day) = mobility_countries[,1]
  
  #Return outputs
  return(list(mobility_countries, sequencing_rate_data, tats_by_day, mobility_data_wide))
}

country_data = getCountryData()
saveRDS("country_data","country_data.rds")



