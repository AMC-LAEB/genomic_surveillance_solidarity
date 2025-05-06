library(tidyr)
library(doParallel)

n_sim = 10000
beta = 0.2
nday = 500
country_data = readRDS("country_data.rds")

# Get matrix P_ij
mobilityMat = country_data[[4]]

# Get country populations (to weight probability for the index locations)
country_populations = country_data[[1]]$Population

# Get origin countries for all simulations
index_countries = sample(1:ncol(mobilityMat),n_sim,replace=TRUE,prob=country_populations)

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


for (r0 in c(1.2,1.3,1.6,2)){
  my.cluster <- parallel::makeCluster(
    50, 
    type = "PSOCK"
  )

  doParallel::registerDoParallel(cl = my.cluster)

  q = foreach(
    i = index_countries
  ) %dopar%{
    simulate_System(i, beta*r0, beta, mobilityMat, nday, 1)
  }
  
  saveRDS(q,paste0("results_",r0,".RDS"))
  
  q = foreach(
    i = index_countries
  ) %dopar%{
    simulate_System(i, beta*r0, beta, mobilityMat*3, nday, 1)
  }
  
  saveRDS(q,paste0("results_",r0,"_fast.RDS"))
  
  q = foreach(
    i = index_countries
  ) %dopar%{
    simulate_System(i, beta*r0, beta, mobilityMat/3, nday, 1)
  }
  
  saveRDS(q,paste0("results_",r0,"_slow.RDS"))
}


