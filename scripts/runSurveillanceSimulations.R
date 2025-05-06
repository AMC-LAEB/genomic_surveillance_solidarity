nday = 500

### Simulate wildtype epidemics by country 
simulateWildtype <- function(beta, N, prevalence, ndays = 500){
  
  gamma = 0.2
  
  out = matrix(0,ncol=6,nrow=ndays)
  out_2 = rep(0,ndays)
  
  state_vec = c(N-N*prevalence, N*prevalence, 0)
  
  for (i in 1:ndays){
    
    new_infections = 0
    
    for (k in 1:10){
      
      S = state_vec[1]
      I = state_vec[2]
      R = state_vec[3]
      
      SI = beta * S * I / N * 0.1
      IR = gamma * I * 0.1
      
      new_infections = new_infections + SI
      
      state_vec[1] = state_vec[1] - SI
      state_vec[2] = state_vec[2] + SI - IR
      state_vec[3] = state_vec[3] + IR
      
    }
    
    state_vec[state_vec < 0] = 0
    out[i,] = state_vec
    out_2[i] = new_infections
    
  }
  
  return(list(out, out_2))
}


### Get sequencing rates by country depending on the strategy
getSequencingRates <- function(min_seqrate, max_seqrate, phi){
  
  seqrates_all =list() # List containing TAT-specific sequencing rates for each country
  
  countrydata = readRDS("country_data.rds") # Dataframe containing info on e.g. sequencing rates for each country
  
  countries = countrydata[[1]][,1]
  
  seq_changed = rep(F,length(countries))
  
  for (ctr in 1:length(countries)){
    
    country = countries[ctr]
    
    #Get country's population size and sequencing rate
    population_size = countrydata[[1]][ctr,]$Population
    sequencing_rate = countrydata[[1]][ctr,]$Seqrate
    
    
    tat_by_day = countrydata[[3]][[ctr]] # Each value corresponds to the proportion of sequences that have TAT in days equal to the index - 1 
    
    seqrate_by_tat = tat_by_day * sequencing_rate # Compute turnaround-time specific sequencing rates
    
    # To account for delays between the sequence readout and submission to GISAID, we assume that the true turnaround time is equal to \phi multiplied
    # by the estimated turnaround time. To account for this, we modify the distribution of turnaround time-specific sequencing rates.
    seqrate_by_tat_new = rep(0,500)
    for (i in 1:nday){
      idxes = pmax(1,floor((1:nday)*phi+0.5))
      seqrate_by_tat_new[idxes[i]] = seqrate_by_tat_new[idxes[i]] + seqrate_by_tat[i]
    }
    seqrate_by_tat = seqrate_by_tat_new
    
    # Get ratio between the country's sequencing rate and the max sequencing rate for this strategy
    ratio = sequencing_rate / max_seqrate
    
    # If this ratio exceeds one, we need to cap that country's sequencing rate. We do so by dividing the vector
    # of turnaround time-specific sequencing rates by the ratio; hence, we evenly reduce sequencing rate across all values 
    # of turnaround time
    
    if (ratio > 1){
      seqrate_by_tat = seqrate_by_tat / ratio
    }
    
    # Now, compute the difference between the minimum sequencing rate and the sequencing rate 
    # Specifically, we look at the sequencing rate with turnaround time <= 14 days, as this is what the minimum sequencing rate requires
    diff = min_seqrate - sum(seqrate_by_tat[1:15]) 
    
    # If there is a deficit, we ensure that the country does attain the minimum sequencing rate by adding the deficit to the 
    # sequencing rate at 1 day's turnaround time (which are collated every fortnight)
    if (diff > 0){
      seq_changed[ctr] = T
      seqrate_by_tat[1] = seqrate_by_tat[1] + diff
    }
    
    # Now, see what the sequencing rate is for turnaround times greater than 14 days
    seqrate_rest = sum(seqrate_by_tat[16:length(seqrate_by_tat)])
    
    # If this sequencing rate is smaller than the deficit (from above) and there is a deficit, we assume that all sequencing output at a greater 
    # turnaround time has had its turnaround time reduced to 14 days. The remaining deficit is added as de novo sequencing capacity
    if (seqrate_rest <= diff & diff > 0){
      seqrate_by_tat[16:length(seqrate_by_tat)] = 0
    }
    
    # If this sequencing rate is greater than the deficit and there is a deficit, then the sequencing rates at greater than 14 day turnaround time
    # are uniformly reduced such that there is no absolute change in sequencing rate (establishing the minimum capacity is done by reducing turnaround time 
    # for the sequencing rates corresponding to a greater turnaround time)
    if (seqrate_rest > diff & diff > 0){
      seqrate_by_tat[16:length(seqrate_by_tat)] = seqrate_by_tat[16:length(seqrate_by_tat)] * (1 - diff/seqrate_rest)
    }
    
    seqrates_all[[ctr]] = seqrate_by_tat #Save turnaround time-specific sequencing rates
  }
  return(list(do.call(rbind,seqrates_all),seq_changed))
}


### Get all wildtype epidemic trajectories
getAllWildtypeTrajectories <- function(prevalence, beta, ndays = 1000){
  
  countrydata = readRDS("country_data.rds") # Dataframe containing info on e.g. sequencing rates for each country
  
  population_sizes = countrydata[[1]]$Population
  
  out = matrix(0, length(population_sizes), ndays)
  
  for (i in 1:length(population_sizes)){
    
    out[i,] = simulateWildtype(beta, population_sizes[i], prevalence, ndays)[[2]]
  }
  
  return(out)
}

runAllSurveillanceSimulations <- function(idx){
  
  library(Rcpp)
  sourceCpp("surveillanceSimulations.cpp")
  
  ### Run genomic surveillance simulations
  
  ## Values of variant Re
  variant_r0s = c(1.2,1.3,1.6,2)
  
  ## Define different sequencing strategies
  min_seqrates = c(0,0,2,2)
  max_seqrates = c(Inf,30,Inf,30)
  
  ## Wildtype incidences at variant introduction
  base_prevs = c(0.001,0.002,0.005,0.02)
  
  ## Wildtype transmission rates
  wildtype_betas = c(0.2,0.21,0.22,0.2)
  
  ## Travel rates
  travel_rates = c("mean","fast","slow")
  
  ## True turnaround time as time factor between collection and submission
  phis = c(0.25,0.5,1)
  
  params = expand.grid(1:4,1:4,1:3,1:3)
  
  min_seqrate = min_seqrates[params[idx,1]]
  max_seqrate = max_seqrates[params[idx,1]]
  phi = phis[params[idx,3]]
  seqrates_all = getSequencingRates(min_seqrate, max_seqrate, phi)
  
  base_prev = base_prevs[params[idx,2]]
  wt_beta = wildtype_betas[params[idx,2]]
  wildtype_epidemics = getAllWildtypeTrajectories(base_prev, wt_beta)
  
  result_matrices =  vector(mode="list",length=length(variant_r0s))
  
  travel_rate_idx = params[idx,4]
  
  countrydata = readRDS("country_data.rds") # Dataframe containing info on e.g. sequencing rates for each country
  population_sizes = countrydata[[1]]$Population
    
  for (r0_idx in 1:length(variant_r0s)){
    
    if (travel_rate_idx == 1){
      path = paste0("simulations/results_",variant_r0s[r0_idx],".RDS")
    }
    if (travel_rate_idx == 2){
      path = paste0("simulations/results_",variant_r0s[r0_idx],"_fast.RDS")
    }
    if (travel_rate_idx == 3){
      path = paste0("simulations/results_",variant_r0s[r0_idx],"_slow.RDS")
    }
    
    simulations = readRDS(path)
    
    results = as.data.frame(getTimeToDetection(simulations,
                                               T,
                                               population_sizes,
                                               seqrates_all[[1]],
                                               wildtype_epidemics,
                                               seqrates_all[[2]],
                                               14))
    
    rm(simulations)
    gc()
    
    colnames(results) = c(
      "detection_day",
      "detection_country",
      "detection_infections",
      "threshold_day",
      "threshold_country",
      "threshold_infections",
      "onset_country")
    
    results$variant_r0 = variant_r0s[r0_idx]
    results$strategy = paste0("min",min_seqrates[params[idx,1]],"max",max_seqrates[params[idx,1]])
    results$base_prev = base_prevs[params[idx,2]]
    results$wildtype_beta = wildtype_betas[params[idx,2]]
    results$phi = phis[params[idx,3]]
    results$travel_rate = travel_rates[params[idx,4]]
    result_matrices[[r0_idx]] = results
  }
  
  detection_df = do.call(rbind,result_matrices)
  saveRDS(detection_df,paste0("detection_outputs/GS_sim_",idx,".rds"))
}


