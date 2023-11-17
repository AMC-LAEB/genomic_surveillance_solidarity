
### Read genomic surveillance outputs
detection_df = readRDS("detection_outputs.rds")
detection_df = detection_df[detection_df$travel_rate=="mean",]


### Only retain main scenarios from main text
detection_df$scenario = NA
detection_df[detection_df$variant_r0==1.2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.001,]$scenario = 1.2
detection_df[detection_df$variant_r0==1.3 & detection_df$wildtype_beta==0.21 & detection_df$base_prev==0.002,]$scenario = 1.3
detection_df[detection_df$variant_r0==1.6 & detection_df$wildtype_beta==0.22 & detection_df$base_prev==0.005,]$scenario = 1.6
detection_df[detection_df$variant_r0==2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.02,]$scenario = 2

detection_df = detection_df[!(is.na(detection_df$scenario)),]

### Get country metadata
cd = readRDS("country_data.rds")[[1]]
continents = cd$Continent
countries = cd$Country
detection_df$onsetcontinent = continents[detection_df$onsetcountry]


#Initialize output
all_ctrydfs = data.frame(country=NA,mean=NA,continent=NA,strat=NA,r0=NA)


idx = 1

for (r0 in c(1.2,1.3,1.6,2)){
  
  results = readRDS(paste0("results_",r0,".RDS")) # Read simulation output
  
  for (i in 1:length(countries)){

    ctry = countries[i]
    
    # Get country-specific onset timings
    onsets = unlist(lapply(1:10000,function(x)which(results[[x]][i,]>0)[1]))
    
    if (length(which(!(is.na(onsets)))) == 0){next}
    
    for (strategy in c("min0maxInf","min2maxInf","min0max30","min2max30","double")){
      
      # Get global detection days
      global_detection = detection_df[detection_df$strategy==strategy & detection_df$variant_r0 == r0 & detection_df$phi==1 & detection_df$travel_rate == "mean",]$detection_day
      diffs = onsets - global_detection # Get time between global detection and local onset
      mn = mean(diffs,na.rm=T) # Compute mean
      all_ctrydfs[idx,] = c(countries[i], mn, continents[i], strategy, r0) # Save 
      idx = idx + 1
    }
  }
}


all_ctrydfs[,2] = as.numeric(all_ctrydfs[,2])

saveRDS(all_ctrydfs,"lead_times.rds")


