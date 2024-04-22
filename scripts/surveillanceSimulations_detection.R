library(tidyr)
library(doParallel)
library(doSNOW)
library(bigmemory)

nsim = 10000
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

binomCI <- function(n, N, p, conf.level){
  pval = pbinom(n - 1, N, p, lower.tail = FALSE)
  CI = c(p.L(n, N, 1 - conf.level), 1)
  return(list(pval, CI))
}


p.L <- function(x, n, alpha) {
  if (x == 0) 
    0
  else qbeta(alpha, x, n - x + 1)
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


### Get sequencing rates by country depending on the strategy
getSequencingRates <- function(min_seqrate, max_seqrate, phi, double){
  
  seqrates_all =list() # List containing TAT-specific sequencing rates for each country
  
  countrydata = readRDS("country_data.rds") # Dataframe containing info on e.g. sequencing rates for each country
  
  countries = countrydata[[1]][,1]
  
  for (ctr in 1:length(countries)){
    
    country = countries[ctr]
    
    #Get country's population size and sequencing rate
    population_size = countrydata[[1]][ctr,]$Population
    sequencing_rate = countrydata[[1]][ctr,]$Seqrate
    
    #For strategy where a country's sequencing rate is doubled
    if (double){sequencing_rate = sequencing_rate * 2}
    
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
    
    # If there is a deficit, we ensure that the country does attain the minimum sequencing rate by adding the deficit to the sequencing rate at 14 day's turnaround time
    if (diff > 0){
      seqrate_by_tat[15] = seqrate_by_tat[15] + diff
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
      seqrate_by_tat[16:length(seqrate_by_tat)] = seqrate_by_tat[16:length(seqrate_by_tat)] * (1- diff/seqrate_rest)
    }
    
    seqrates_all[[ctr]] = seqrate_by_tat #Save turnaround time-specific sequencing rates
  }
  return(seqrates_all)
}

getTimesSim <- function(out, nsim, wt_beta, base_prev, min_seqrate, max_seqrate, phi, double, isThreshold){
  
  nday = 500
  
  countrydata = readRDS("country_data.rds") #Get metadata
  countries = countrydata[[1]][,1]
  
  ### Simulate wildtype epidemics for all countries
  wildtype_epidemics = getAllWildtypeTrajectories(base_prev, wt_beta)
  
  ### Get turnaround time-specific sequencing rates for the strategy
  seqrates_all = getSequencingRates(min_seqrate, max_seqrate, phi, double)
  
  getStatsSingle <- function(simidx, isThreshold){
    
    isDetected = F
    
    day_of_minimum_threshold = NA
    country_of_minimum_threshold = NA
    incidence_on_day_of_minimum_threshold = NA
    
    results = out[(((simidx-1)*196)+1):((simidx*196)),] # Extract this simulation from matrix of simulations
    
    onsets = apply(results,1,function(x)which(x>0)[1]) # Get country-specific onsets for this simulation
    
    index_country = which.min(onsets) # Get simulation index country
    
    mt_count_matrix = matrix(0,196,nday) # Keep track of number of variant sequences
    wt_count_matrix = matrix(0,196,nday) # Keep track of number of wildtype sequences
    
    p_value_matrix = matrix(1,196,nday) # Keep track of p-values for binomial tests that the variant proportion exceeds a threshold
    
    for (day in 1:nday){ # Iterate through all days
      
      for (ctr in which(onsets<=day)){ # Only look at countries that have seen the variant arrive already
        
        onset_day = which(results[ctr,]>0)[1] # Get country's onset day
        
        plus_10_day = which(results[ctr,]>=10)[1] # Get first day variant incidence exceeds or is equal to 10
        
        seqrate_by_tat = seqrates_all[[ctr]] # Get sequencing rates by TAT
        
        population_size = countrydata[[1]][ctr,]$Population # Country's population size
        
        prop_thresh = 0.001 * (1e8 / population_size) # Get threshold that the variant proportion should exceed for this country
        
        if (prop_thresh > 0.01){prop_thresh = 0.01}
        
        wildtype_epidemic = wildtype_epidemics[ctr,] # Get country's wildtype epidemic
        
        ndays_sim = nday - plus_10_day + 1 # Get number of days since the first day incidence is at or above 10
        
        variant_proportion = rep(0,nday) # Initialize variant proportion
        
        base_inc = wildtype_epidemic[1] # Get initial wildtype incidence
        
        wildtype_cases = rep(0,nday) # Initialize wildtype epidemic for this specific simulation
        
        if (!(is.na(plus_10_day))){
          
          # Until the first day variant incidence is at or above 10, wildtype incidence is equal to initial wildtype incidence
          wildtype_cases[1:(plus_10_day-1)] = base_inc 
          # From the first day variant incidence is at or above 10, wildtype incidence is given by the simulated wildtype dynamics
          wildtype_cases[plus_10_day:nday] = wildtype_epidemic[1:ndays_sim]
          # Compute the variant proportion given the simulated variant and wildtype epidemic dynamics
          variant_proportion = results[ctr,]/(results[ctr,] + wildtype_cases)
          
        } else {
          
          # If variant incidence never reached 10, wildtype incidence is always equal to the initial wildtype incidence
          wildtype_cases[] = base_inc
          # Compute the variant proportion given the simulated variant and wildtype epidemic dynamics
          variant_proportion = results[ctr, ]/(results[ctr, ]+ wildtype_cases)
        }
        
        # Get time between the current day and the day of the first case
        time_since_onset = day - onset_day + 1
        
        # Convert from S/M/wk to S/d
        true_seqrates = seqrate_by_tat[1:time_since_onset]*(population_size/1e6)/7
        
        # For all days between the onset day and the current day, get the Poisson-distributed number of sequences that will be available on the current day
        # Because this array is reversed, the first element corresponds to the number of sequences  that become available today and were collected on the first day since the variant's onset
        # The last element corresponds to the number of sequences that became available today for which samples were collected today 
        number_of_sequences = rpois(time_since_onset,rev(true_seqrates))
        # Cap number of sequences at true number of cases
        number_of_sequences = pmin(number_of_sequences, wildtype_cases[onset_day:day])
        # Using these number of sequences, simulate the binomial number of variant sequences whose results became available on the current day, for each day in the past since the variant's onset
        nseq_mt_on_day = rbinom(time_since_onset, number_of_sequences, variant_proportion[onset_day:day])
        nseq_wt_on_day = number_of_sequences - nseq_mt_on_day
        
        # Save the number of variant and wildtype sequences returned on each day
        mt_count_matrix[ctr,onset_day:day] = mt_count_matrix[ctr,onset_day:day] + nseq_mt_on_day
        wt_count_matrix[ctr,onset_day:day] = wt_count_matrix[ctr,onset_day:day] + nseq_wt_on_day
        
        #If this is the index day in the week, go through all previous weeks to see if the variant proportion has exceeded 0.1% equiv at any previous week
        
        if (day %% 7 == 0){
          all_days = which(1:day%%7==0) # Get first days of weeks
          for (test_day in all_days){ # Check for all first days of weeks in the past
            mt_tot = sum(mt_count_matrix[ctr,(test_day-6):test_day],na.rm=T) # Get number of variant sequences collected in that week
            wt_tot = sum(wt_count_matrix[ctr,(test_day-6):test_day],na.rm=T) # Get number of wildtype sequences collected in that week
            
            if (mt_tot < 1){next}
            
            
            if (!(isThreshold)){
              p_value = binomCI(mt_tot, mt_tot+wt_tot, prop_thresh, 0.05)[[1]] # Get p-value of one-sided binomial test that variant prop exceeds thresh
              p_value_matrix[ctr,day] = p_value # Save p-value
              
              if (p_value < 0.05){
                day_of_minimum_threshold = day # Save day estimated variant proportion exceeds threshold
                country_of_minimum_threshold = ctr # Save country estimated variant proportion exceeds threshold
                incidence_on_day_of_minimum_threshold = sum(results[,1:day],na.rm=T) # Total number of global cases on day variant proportion exceeds thresh
                isThreshold = T
              }
            }
          }
        }
      }
      
      total_seqs = sum(mt_count_matrix,na.rm=T) # Get total number of variant sequences on this day
      
      #Save the day of detection, incidence on that day, and the country
      if (total_seqs >= 1 & !(isDetected)){
        day_of_detection = day # Save day of detection
        country_of_first_detection = which.max(rowSums(mt_count_matrix)) # Save country that detected first
        incidence_on_day_of_detection = sum(results[,1:day],na.rm=T) # Total number of global cases on day of detection
        isDetected = T
      }
      
      if (isDetected & isThreshold){ # Return if has been detected and seen to exceed threshold in at least one country
        return(c(day_of_detection,
                 country_of_first_detection,
                 incidence_on_day_of_detection,
                 day_of_minimum_threshold,
                 country_of_minimum_threshold,
                 incidence_on_day_of_minimum_threshold,
                 index_country))
      }
    }
    if (isDetected){ 
      return(c(day_of_detection,
               country_of_first_detection,
               incidence_on_day_of_detection,
               NA,
               NA,
               NA,
               index_country))
    }
    return(c(NA,NA,NA,NA,NA,NA,index_country))
  }
  
  #Iterate through all simulations
  detection_outputs = matrix(NA,nsim,7)
  for (i in 1:nsim){
    detection_outputs[i,] = getStatsSingle(i, isThreshold)
  }
  
  return(detection_outputs)
}

### Run genomic surveillance simulations

## Values of variant Re
variant_r0s = c(1.2,1.3,1.6,2)

## Define different sequencing strategies
min_seqrates = c(0,0,2,2,0)
max_seqrates = c(Inf,30,Inf,30,Inf)
double = c(F,F,F,F,T)

## Wildtype incidences at variant introduction
base_prev = c(0.001,0.002,0.005,0.02)

## Wildtype transmission rates
wildtype_betas = c(0.2,0.21,0.22,0.2)

## Travel rates
travel_rates = c("mean","fast","slow")


## True turnaround time as time factor between collection and submission
phis = c(0.25,0.5,1)


### Initialize cluster
my.cluster <- parallel::makeCluster(
  40
)

registerDoSNOW(my.cluster)

### Initialize output
all_results = list() 
result_idx = 1

params = expand.grid(1:5,1:4,1:3)


### Get epidemic simulation output memory-efficient

result_matrices =  vector(mode="list",length=length(variant_r0s))


### Simulate genomic surveillance

for (r0_idx in 1:length(variant_r0s)){
  
  print(paste0("Simulating surveillance for variant R0 ", variant_r0s[r0_idx]))
  
  for (travel_rate_idx in 1:length(travel_rates)){
    
    print(paste0("Simulating surveillance for travel rate ", travel_rates[travel_rate_idx]))
    
    if (travel_rate_idx == 1){
      path = paste0("results_",variant_r0s[r0_idx],".RDS")
    }
    if (travel_rate_idx == 2){
      path = paste0("results_",variant_r0s[r0_idx],"_fast.RDS")
    }
    if (travel_rate_idx == 3){
      path = paste0("results_",variant_r0s[r0_idx],"_slow.RDS")
    }
    
    results = readRDS(path)
    mat = as.big.matrix(do.call(rbind,results))
    mdesc = describe(mat)
    
    GS_results = foreach(
      i = 1:nrow(params)
    ) %dopar%{
      #) %do% {
      require(bigmemory)
      getTimesSim(attach.big.matrix(mdesc),
                  nsim,
                  wildtype_betas[params[i,2]],
                  base_prev[params[i,2]],
                  min_seqrates[params[i,1]],
                  max_seqrates[params[i,1]],
                  phis[params[i,3]],
                  double[params[i,1]],
                  T)
    }
    
    GS_results_out = GS_results
    
    for (i in 1:length(GS_results)){
      
      df = as.data.frame(GS_results[[i]])
      colnames(df) = c(
        "detection_day",
        "detection_country",
        "detection_infections",
        "threshold_day",
        "threshold_country",
        "threshold_infections",
        "onset_country")
      
      df$variant_r0 =variant_r0s[r0_idx]
      if (!(double[params[i,1]])){
        df$strategy = paste0("min",min_seqrates[params[i,1]],"max",max_seqrates[params[i,1]])
      } else {
        df$strategy = "double"
      }
      df$base_prev = base_prev[params[i,2]]
      df$wildtype_beta = wildtype_betas[params[i,2]]
      df$phi = phis[params[i,3]]
      df$travel_rate = travel_rates[travel_rate_idx]
      GS_results_out[[i]] = df
    }
    
    GS_results_df = do.call(rbind,GS_results_out)
    all_results[[result_idx]] = GS_results_df
    result_idx = result_idx + 1
  }
}

detection_df = do.call(rbind,all_results)
saveRDS(detection_df,'detection_outputs_threshold_F.rds')
