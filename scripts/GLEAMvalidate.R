library(tidyr)
library(doParallel)
library(ggplot2)
library(lemon)


### Run epidemic simulations 


# Simulation params
n_sim = 10000
beta = 0.2
nday = 1000

# Country info
country_data = readRDS("country_data.rds")
cd = country_data[[1]]


# Get matrix P_ij
mobilityMat = country_data[[4]]

# Function to simulate epidemic dynamics
simulate_System <- function(indexcountry, alpha, beta, mat, nday, tau){
  
  ncountry = dim(mat)[1]
  
  # Save number of new infections on each day
  infections = matrix(0, ncountry, nday)
  
  # Get population sizes
  n = country_data[[1]]$Population
  
  while (sum(infections,na.rm=T) < 1000){ # Repeat if stochastic extinction occurred 
    
    infections = matrix(0, ncountry, nday)
    
    j = rep(0, ncountry)
    
    j[indexcountry] = 10
    
    s = n
    
    for (i in 1:(nday)){

      if (sum(j,na.rm=T) == 0){break} # Stop if no-one is infectious
      
      for (rep in 1:(1/tau)){
      
        new_infections = rpois(ncountry, tau * alpha * s * j / n)
      
        recoveries = rpois(ncountry, tau * beta * j) 
      
        j_mobility = sign(rowSums(t(t(mat)*j) - (t(mat)*j),na.rm=T)) * rpois(ncountry, tau * (abs(rowSums(t(t(mat)*j) - (t(mat)*j),na.rm=T))))
        s_mobility = sign(rowSums(t(t(mat)*s) - (t(mat)*s),na.rm=T)) * rpois(ncountry, tau * (abs(rowSums(t(t(mat)*s) - (t(mat)*s),na.rm=T))))
      
        infections[,i] = infections[,i] + new_infections
       
        j = j + new_infections - recoveries + j_mobility
      
        s = s - new_infections + s_mobility
      
        s[s<0] = 0
        j[j<0] = 0
      }
    }
  }
  
  mode(infections) = "integer" # Reduce output size 
  
  return(infections)
}

country_name_list = list()
country_name_list[[1]] = c("Cape Verde", "Cabo Verde")
country_name_list[[2]] = c("Côte d'Ivoire", "Cote d'Ivoire","Ivory Coast")
country_name_list[[3]] = c("Czech Republic","Czechia")
country_name_list[[4]] = c("Congo, Dem. Rep.", "Democratic Republic of the Congo","DR Congo","Congo (Democratic Republic of the)")
country_name_list[[5]] = c("Egypt, Arab Rep.", "Egypt")
country_name_list[[6]] = c("Swaziland", "Eswatini")
country_name_list[[7]] = c("Lao PDR", "Laos","Lao People's Democratic Republic")
country_name_list[[8]] = c("Macedonia", "North Macedonia","Macedonia (the former Yugoslav Republic of)")
country_name_list[[9]] = c("Congo, Rep.", "Republic of the Congo","Congo")
country_name_list[[10]] = c("Russian Federation", "Russia")
country_name_list[[11]] = c("St. Maarten", "Sint Maarten")
country_name_list[[12]] = c("Slovak Republic", "Slobakia")
country_name_list[[13]] = c("Korea, Rep.", "South Korea","Korea (Republic of)")
country_name_list[[14]] = c("Virgin Islands (U.S.)", "U.S. Virgin Islands","United States Virgin Islands")
country_name_list[[15]] = c("St. Lucia", "Saint Lucia")
country_name_list[[16]] = c("St. Pierre and Miquelon", "Saint Pierre and Miquelon")
country_name_list[[17]] = c("St. Vincent and Grenadines", "Saint Vincent and Grenadines", "Saint Vincent and the Grenadines")
country_name_list[[18]] = c("St. Kitts and Nevis", "Saint Kitts and Nevis")
country_name_list[[19]] = c("Slovak Republic", "Slovakia")
country_name_list[[20]] = c("Micronesia, Fed. States of", "Micronesia","Micronesia (Federated States of)")
country_name_list[[21]] = c("Korea, Dem. Rep.", "North Korea","Korea (Democratic People's Republic of)")
country_name_list[[22]] = c("Wallis and Futuna Islands", "Wallis and Futuna")
country_name_list[[23]] = c("Curaçao", "Curacao")
country_name_list[[24]] = c("Faeroe Islands", "Faroe Islands")
country_name_list[[25]] = c("São Tomé and Principe", "Sao Tome and Principe")
country_name_list[[26]] = c("St-Barthélemy", "Saint Barthelemy")
country_name_list[[27]] = c("Trinidad and Tobago", "Trinidad & Tobago")
country_name_list[[28]] = c("United States of America", "United States", "USA")
country_name_list[[29]] = c("Kyrgyzstan", "Kyrgyz Republic")
country_name_list[[30]] = c("Iran", "Iran, Islamic Rep.","Iran (Islamic Republic of)")
country_name_list[[31]] = c("Gambia", "Gambia, The")
country_name_list[[32]] = c("Turkey", "Turkiye")
country_name_list[[33]] = c("Venezuela", "Venezuela, RB", "Venezuela (Bolivarian Republic of)")
country_name_list[[34]] = c("Syria", "Syrian Arab Republic")
country_name_list[[35]] = c("Yemen", "Yemen, Rep.")
country_name_list[[36]] = c("Bolivia", "Bolivia (Plurinational State of)")
country_name_list[[37]] = c("Brunei", "Brunei Darussalam")
country_name_list[[38]] = c("United Kingdom", "United Kingdom of Great Britain and Northern Ireland")
country_name_list[[39]] = c("British Virgin Islands", "Virgin Islands (British)")
country_name_list[[40]] = c("Moldova", "Moldova (Republic of)")
country_name_list[[41]] = c("Tanzania", "Tanzania, United Republic of")
country_name_list[[42]] = c("Vietnam", "Viet Nam")
country_name_list[[43]] = c("Bahamas","The Bahamas")
country_name_list[[44]] = c("Palestine","Palestine, State of")
country_name_list[[45]] = c("Macao","Macau")


# Get country names in GLEAM and match to my country names
ctrnames = read.csv("GLEAM_validation/country_names.tsv",sep='\t')
for (i in 1:length(country_name_list)){
  whch_s = which(ctrnames$Name %in% country_name_list[[i]])
  if (length(whch_s)>0){
    ctrnames[whch_s,]$Name = country_name_list[[i]][1]
  }
}

# Run simulations for all the countries in the sensitivity analysis

results = data.frame(index_country=NA,country=NA,gleam=NA,metapop=NA)
rowidx = 1

countries = c("Mali","Oman","Nepal","Jamaica","Uzbekistan","Malaysia","Ecuador","Nicaragua","France","Cameroon")

for (country_idx in 1:length(countries)){
  
  print(paste0("Simulating for ",countries[country_idx]))
  sims_metapop = lapply(1:10,function(x)simulate_System(match(countries[country_idx],colnames(mobilityMat)), 0.2*1.6, 0.2, mobilityMat, 250, 1))
  
  for (idx in 1:nrow(ctrnames)){
    
    ctr = ctrnames[idx,2]
    ctry = ctrnames[idx,]$Name
    cuminc_gleam = read.csv(paste0("GLEAM_validation/",countries[country_idx],"/countries/",ctr,"-0.tsv"),sep='\t')
    day_gleam = which(cuminc_gleam$Cumulative.Median>0.01)[1]
    
    ctr = match(ctry,colnames(mobilityMat))
    cuminc_metapop = do.call(rbind,lapply(1:10,function(x)cumsum(sims_metapop[[x]][ctr,])))
    day_metapop = which(apply(cuminc_metapop,2,median)>(0.01/1000)*cd$Population[ctr])[1]
    
    if (!(is.na(val) & !(is.na(day_gleam)))){
      results[rowidx,] = c(countries[country_idx],ctry,day_gleam,day_metapop)
      rowidx = rowidx + 1
    }
  }
}


results[,3] = as.numeric(results[,3])
results[,4] = as.numeric(results[,4])

saveRDS(results,"GLEAM_validation.rds")


