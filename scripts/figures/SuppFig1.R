simulateEpidemic <- function(state_vec, params, ndays,N, introductionday){
  
  beta_1 = params[1]
  beta_2 = params[2]
  gamma = params[3]
  
  out = matrix(0,ncol=6,nrow=ndays)
  out_2 = matrix(0,ncol=2,nrow=ndays)
  
  for (i in 1:ndays){
    
    if (i == introductionday){
      state_vec[5] = 1
    }
    
    
    S1 = state_vec[1]
    I1 = state_vec[2]
    R1 = state_vec[3]
    S2 = state_vec[4]
    I2 = state_vec[5]
    R2 = state_vec[6]
    
    SI_1 = beta_1 * S1 * I1 / N# * 0.1
    SI_2 = beta_2 * S2 * I2 / N# * 0.1
    IR_1 = gamma * I1 #* 0.1
    IR_2 = gamma * I2 #* 0.1
    
    
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
    
    
    state_vec[state_vec < 0] = 0
    out[i,] = state_vec
    out_2[i,] = c(SI_1,SI_2)
    
  }
  
  return(list(out, out_2))
}

getPredictedTime <- function(f0,s,q,n){
  return((log((q^(-s/n)-1)/f0+1))/s)
}


sequencing_rates = 10^seq(-1,3,0.1)

modeled_time_to_detection = c()
simulated_time_to_detection = c()

ndays = 1000
popsize = 100e6


initial_weekly_wildtype_incidences = c(0.001,0.002,0.005,0.02)
variant_r0s = c(1.2, 1.3, 1.6, 2)
wildtype_betas = c(0.2,0.21,0.22,0.2)

time_df = data.frame(variant_r0=NA,scenario=NA,seqrate=NA,simulated=NA,equation=NA)

for (variant_r0_idx in 1:4){
  for (scenario in 1:4){
    
    variant_beta = variant_r0s[variant_r0_idx]*0.2
    wildtype_beta = wildtype_betas[scenario]
    
    init_I = initial_weekly_wildtype_incidences[scenario] * popsize
    
    state_vec = c(popsize-init_I,init_I,0,100e6,0,0)
    
    epi_params = c(wildtype_beta, variant_beta, 0.2)
    
    sim = simulateEpidemic(state_vec,epi_params,1000,popsize,1)

    incidence_mat = matrix(c(sim[[2]][,1],sim[[2]][,2]),nrow=1000)
    
    for (sequencing_rate in sequencing_rates){
      
      detection_day = median(replicate(100,which(cumprod(dbinom(0,rpois(1000,sequencing_rate*popsize/1e6/7),sim[[2]][,2]/rowSums(sim[[2]][,])))<0.05)[1]),na.rm=T)
      modeled_time_to_detection = detection_day
      simulated_time_to_detection = getPredictedTime((incidence_mat[1,2])/sum(incidence_mat[1,]),variant_beta-wildtype_beta,0.05,popsize/1e6*sequencing_rate/7)
      time_df = rbind(time_df,c(variant_r0s[variant_r0_idx],scenario,sequencing_rate,modeled_time_to_detection,simulated_time_to_detection))
    }
  }
}



thm = theme(axis.text = element_text(color="black",size=5),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position='none',
            axis.title= element_text(size=5),
            strip.background=element_rect(colour="white",fill="white"),
            strip.text = element_text(size = 5),
            legend.text = element_text(size = 5),
            legend.title=element_text(size=5),
            legend.spacing.x = unit(0.05,"cm"),
            legend.key=element_rect(fill="white"))


time_df$scenario = factor(time_df$scenario, 
                                       levels=c(1,2,3,4),
                                       labels=c(expression(atop(paste("wt ",italic(R[e])," = 1"),paste("wt prevalence = 0.1%"))),
                                                expression(atop(paste("wt ",italic(R[e])," = 1.05"),paste("wt prevalence = 0.2%"))),
                                                expression(atop(paste("wt ",italic(R[e])," = 1.1"),paste("wt prevalence = 0.5%"))),
                                                expression(atop(paste("wt ",italic(R[e])," = 1"),paste("wt prevalence = 2%")))))

time_df$variant_r0 = factor(time_df$variant_r0, 
                          levels=c(1.2,1.3,1.6,2),
                          labels=c(expression(paste("Variant ",italic(R[e])," = 1.2")),
                                   expression(paste("Variant ",italic(R[e])," = 1.3")),
                                   expression(paste("Variant ",italic(R[e])," = 1.6")),
                                   expression(paste("Variant ",italic(R[e])," = 2"))))




time_df = time_df[!(is.na(time_df$variant_r0)),]
ggplot(time_df,aes(x=simulated,y=equation)) + geom_point(cex=0.5,col='red') + 
  facet_rep_grid(scenario~variant_r0,scales='free',labeller=label_parsed) + 
  geom_abline(linewidth=0.5,lty=2,col='grey') + thm + xlab("Detection day (simulations)") +  ylab("Detection day (equation)")

ggsave("figures/SuppFig1.pdf",width=120,units="mm",height=100)
