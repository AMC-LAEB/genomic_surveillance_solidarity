library(ggplot2)
library(npreg)
library(patchwork)
library(latex2exp)
library(scales)

# Read country sequencing rate information

countrydata = readRDS("country_data.rds")
cd = countrydata[[2]]

high_median = median(cd[cd$Income=="High income",]$Seqrate,na.rm=T)
low_median = median(cd[cd$Income=="Low income",]$Seqrate,na.rm=T)


### Simulates the SIR epidemic dynamics of two viruses

simulateEpidemic <- function(state_vec, params, ndays, N){
  
  beta_1 = params[1]
  beta_2 = params[2]
  gamma = params[3]
  
  out = matrix(0,ncol=6,nrow=ndays)
  out_2 = matrix(0,ncol=2,nrow=ndays)
  
  for (i in 1:ndays){
    
    variant_infections = 0
    wildtype_infections = 0
    
    for (k in 1:10){
      
      S1 = state_vec[1]
      I1 = state_vec[2]
      R1 = state_vec[3]
      S2 = state_vec[4]
      I2 = state_vec[5]
      R2 = state_vec[6]
      
      SI_1 = beta_1 * S1 * I1 / N * 0.1
      SI_2 = beta_2 * S2 * I2 / N * 0.1
      IR_1 = gamma * I1 * 0.1
      IR_2 = gamma * I2 * 0.1
      
      
      dS1 <- -SI_1
      dI1 <- SI_1 - IR_1
      dR1 <- IR_1
      dS2 <- -SI_2
      dI2 <- SI_2 - IR_2
      dR2 <- IR_2
      
      state_vec[1] = state_vec[1] + dS1
      state_vec[2] = state_vec[2] + dI1
      state_vec[3] = state_vec[3] + dR1
      state_vec[4] = state_vec[4] + dS2
      state_vec[5] = state_vec[5] + dI2
      state_vec[6] = state_vec[6] + dR2
      
      wildtype_infections = wildtype_infections + SI_1
      variant_infections = variant_infections + SI_2
      
    }
    
    state_vec[state_vec < 0] = 0
    out[i,] = state_vec
    out_2[i,] = c(wildtype_infections, variant_infections)
    
  }
  
  return(list(out, out_2))
}

# Given a sequencing rate in units of S/d and a matrix of wildtype and variant infections through time, compute
# the day by which the variant will have been detected with 95% confidence

getDayOfDetection <- function(seqRatePerDay, incidence_matrix){
  nseq_per_day = floor(pmin(rpois(1000,seqRatePerDay),rowSums(incidence_matrix))) # Sample Poisson-valued number of samples to sequence, capped at true number of infections
  variant_prop = incidence_matrix[,2]/rowSums(incidence_matrix) # Compute variant proportion through time
  day_of_detection = which(cumprod(dbinom(0,nseq_per_day,variant_prop))<0.05)[1] 
  return(day_of_detection)
}


## Define simulation parameters

gamma = 0.2 # Recovery rate
popsize= 1e8 # Population size
tat = 14 # Turnaround time 

# Sequencing rates
sequencing_rates_smwk = sort(c(10^seq(-2,3,0.02),low_median,high_median)) 


# Simulation parameters
variant_r0s = c(1.2,1.3,1.6,2) # Variant R_e
wildtype_r0s = c(1,1.05,1.1,1) # Wildtype R_e
wildtype_prev = c(0.001,0.002,0.005,0.02) # Initial wildtype prevalence


all_times_to_detection = list()
all_effective_increases = list()

## Simulate time to variant detection

idx = 1

for (i in 1:length(variant_r0s)){ # Iterate through variant R_e 
  
  print(paste0("Simulating for variant Re ",variant_r0s[i]))
  
  for (j in 1:length(wildtype_r0s)){ # Iterate through scenario of variant emergence (wt initial Re and prevalence)
     
    init_I = wildtype_prev[j] * popsize # Compute initial number of wildtype infecteds
    
    state_vec = c(popsize-init_I, init_I, 0, popsize, 1, 0) # Initialize SIR state vector (SIR[wildtype],SIR[variant])
    
    epi_params = c(wildtype_r0s[j]*gamma, variant_r0s[i]*gamma, gamma) # Initialize SIR parameter vector (wildtype beta, variant beta, gamma)
    
    incidence = simulateEpidemic(state_vec, epi_params, 1000, popsize)[[2]] # Simulate epidemics and extract wildtype and variant incidence
    
    # For each sequencing rate, simulate Poisson-valued number of sequences per day, and, given the computed variant
    # proportion through time, compute the day on which the probability that the variant has not been detected declines below 0.05
    simulations = replicate(100,sapply(sequencing_rates_smwk,function(x)getDayOfDetection(x*popsize/1e6/7,incidence)))
    predicted_time_to_detection = apply(simulations,1,median,na.rm=T)
    
    # Fit a smoothing spline to the relationship between sequencing rate and time to detection
    time_fit = ss(log(sequencing_rates_smwk,10),  predicted_time_to_detection,all.knots = T,lambda = 1e-5) 
    
    # Fit a smoothing spline to the relationship between sequencing rate and cumulative infections by time of detection
    inc_fit = ss(1:length(incidence[,2]),log(cumsum(incidence[,2])),all.knots=T,lambda = 1e-10) 
    
    # Get time of detection when including TAT
    predicted_time_to_detection_tat = predict.ss(time_fit,x=log(sequencing_rates_smwk,10))$y+tat 
    # Get time of detection when including TAT if sequencing rate was increased by 1 S/M/wk
    predicted_time_to_detection_tat_delta = predict.ss(time_fit,x=log(sequencing_rates_smwk+1,10))$y+tat 
    #For each sequencing rate, compute the change in time to variant detection if the sequencing rate were increased by 1 S/M/wk
    predicted_time_to_detection_dn = predicted_time_to_detection_tat - predicted_time_to_detection_tat_delta
    
    # Get incidence at detection when including TAT
    predicted_inc_at_detection = exp(predict.ss(inc_fit,x=predicted_time_to_detection_tat)$y) 
    # Get reduction in incidence at detection when increasing by 1 S/M/wk
    predicted_inc_at_detection_dn = predicted_inc_at_detection - exp(predict.ss(inc_fit,x=predicted_time_to_detection_tat_delta)$y) 
    
    
    ### Compute equivalent increases in turnaround time
    
    detection_days = predicted_time_to_detection[sequencing_rates_smwk > 1 & sequencing_rates_smwk < 100] 
    seqrates_log = log(sequencing_rates_smwk[sequencing_rates_smwk > 1 & sequencing_rates_smwk < 100],10)
    mod = lm(seqrates_log~detection_days)
    equivalent_increases = 10^predict.lm(mod,new=data.frame(detection_days=c(1)))/10^predict.lm(mod,new=data.frame(detection_days=c(1:30)))
    
    ### Collate output
    
    all_effective_increases[[idx]] = data.frame(variant_r0=variant_r0s[i],
                                                wildtype_r0=wildtype_r0s[j],
                                                equivalent_increase = equivalent_increases, 
                                                tat_reduction=1:30, 
                                                initial_wildtype_prevalence=wildtype_prev[j])
    

    
    virus_df = data.frame(variant_r0=variant_r0s[i], 
                          wildtype_r0=wildtype_r0s[j], 
                          seqrate=sequencing_rates_smwk, 
                          initial_wildtype_prevalence=wildtype_prev[j],  
                          day_of_detection = predicted_time_to_detection_tat, 
                          day_of_detection_dn = predicted_time_to_detection_dn, 
                          inc_at_detection = predicted_inc_at_detection, 
                          inc_at_detection_dn = predicted_inc_at_detection_dn)
    
    all_times_to_detection[[idx]] = virus_df
    
    idx = idx + 1
    
  }
}

### Join all outputs
time_to_detection_df = do.call(rbind,all_times_to_detection)
time_to_detection_df$variant_r0 = factor(time_to_detection_df$variant_r0,levels=variant_r0s)
time_to_detection_df = time_to_detection_df[!(is.na(time_to_detection_df$variant_r0)),]
effective_increase_df = do.call(rbind, all_effective_increases)

saveRDS(time_to_detection_df, "time_to_detection.rds")
saveRDS(effective_increase_df, "effective_increase.rds")

