library(ggnewscale)
library(viridis)
library(ggplot2)
library(patchwork)
library(scales)

countrydata = readRDS("../outputs/country_data.rds")
population_sizes = countrydata[[1]]$Population
ctrs = countrydata[[1]]$Country

# This function returns the travel rate from each country in the input vector to the most connected country
# that has sequencing rate >= 10 S/M/wk
getMaxTravelRate <- function(ctries){
  popsizes = countrydata[[1]]$Population[ctries]
  ctries = ctrs[ctries]
  f = read.csv("../outputs/GLEAM_validation/GTM_data.tsv",sep='\t')
  u = reshape2::dcast(f[f$year==2016,][,c(1,2,6)],source_name~target_name)[,2:197]
  ctrs_gtm = colnames(u)
  ctrs_hi_seqrate = ctrs[which(countrydata[[1]]$Seqrate>10)]
  hi_seqrate_gtm = match(ctrs_hi_seqrate,ctrs_gtm)
  rates = sapply(match(ctries,ctrs_gtm),function(x)max(u[hi_seqrate_gtm,x],na.rm=T))
  return(rates/popsizes/365)
}

# Get the countries that have a sequencing rate <0.1 S/M/wk and population >= 1M
low_ctrs =sort(getMaxTravelRate(which(countrydata[[1]]$Seqrate<0.1 & population_sizes>1e6)))


# Get population and epidemiological parameters for epidemic simulation
popsize = 100e6

initial_weekly_wildtype_incidences = c(0.001,0.002,0.005,0.02)
variant_r0s = c(1.2, 1.3, 1.6, 2)
wildtype_betas = c(0.2,0.21,0.22,0.2)


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


# This function returns a list containing 1) a figure showing the expected time to detection based on local vs sentinel detection,
# a dataframe detailing the time to detection at sentinel hubs, and a dataframe detailing time to detection through local sequencing
plotSeqrateVsTravel <- function(variant_beta,wildtype_beta,init_I,q){
  
  state_vec = c(popsize-init_I,init_I,0,popsize,1,0) # Initialize the vector of mt and wt infections
  epi_params = c(wildtype_beta, variant_beta, 0.2)
  sim = simulateEpidemic(state_vec,epi_params,1000,popsize) # Simulate epidemic
  incidence_mat = matrix(c(sim[[2]][,1],sim[[2]][,2]),nrow=1000) # Get prevalence matrix
  df = data.frame(travel=NA,psamp=NA,day=NA,country=NA,variant_beta=NA)
  travelprops = 10^seq(-7,-1,0.1) # Range of daily outward travel rates 
  # Loop through travel rates and ascertainment rates
  for (i in c(low_ctrs,travelprops)){
    for (j in c(1,0.1,0.01)){
      travel_time_to_d = getExportDate(i*j,round(sim[[1]][,5]),q) #Get predicted detection time at sentinel hub 
      df = rbind(df,c(i,j,travel_time_to_d,i%in%low_ctrs,variant_beta))
    }
  }
  
  # Now, get expected detection times based on local sequencing. Because of 14-day sequencing regularity, the expected turnaround time for detection is 7 days
  seqrates = c(0.1,0.2,0.5,1)
  seq_time_to_d = sapply(seqrates,function(x)median(replicate(100,which(cumprod(dbinom(0,rpois(1000,x*popsize/1e6/7),sim[[2]][,2]/rowSums(sim[[2]][,])))<q)[1]))+7)
  seq_df = data.frame(seqrate=seqrates,detection=seq_time_to_d,variant_beta=variant_beta)
  
  return(list(ggplot(df[complete.cases(df),],aes(x=travel,y=day,group=factor(psamp),color=factor(psamp))) + geom_line(lty=1,lwd=.4) +
           theme_bw() + xlab("Daily outward travel rate") + ylab("Median day of detection\nin sentinel country") + 
           scale_color_manual(values=c('darkblue','darkgreen','#5EDC1F'),name='Ascertainment rate',guide=guide_legend(order = 1, direction = "horizontal")) + 
           geom_rug(data=data.frame(rate=low_ctrs),aes(x=rate),inherit.aes=F,size=.3) + 
           scale_x_log10(limits=10^c(-7,-3),labels=trans_format("log10", math_format(10^.x))) + 
           new_scale_colour() + 
           geom_hline(data=seq_df,aes(yintercept=seq_time_to_d,color=factor(seqrate)), linetype=2,linewidth=0.4)+ #xlim(10^c(-7,-1)) +
           theme(axis.text=element_text(color='black'),
                 axis.title=element_text(size=8),
                 legend.direction='vertical',
                 legend.box = "vertical",
                 legend.text = element_text(size=8),
                 legend.spacing.y = unit(-.1, "cm"),
                 legend.key.size = unit(0.3,"cm"),
                 legend.title= element_text(size=8),
                 legend.position='top') + 
           scale_color_manual(values= colors <-magma(6)[2:5],guide=guide_legend(order = 2, direction = "horizontal"),name = "Sequencing rate (S/M/wk)"),
           df,
           seq_df))

}

# This function implements the mathematical model to estimate the expected time until detection at a travel hub,
# given n (the ascertainment rate x the outward travel rate from the origin to sentinel country), prev (the prevalence vector in the origin country),
# and q, the confidence level
getExportDate <- function(n,prev,q){
  mod = lm(log(prev[20:100])~c(20:100))
  c = exp(mod$coefficients[1])
  r = mod$coefficients[2]
  return(as.vector((1/r) * log(1-((log(q)*(exp(r)-1))/(n*c)))))
}

plots = list()
dfs_seq = list()
dfs_travel = list()
idx = 1

for (variant_r0_idx in 1:4){
    scenario = variant_r0_idx
    variant_beta = variant_r0s[variant_r0_idx]*0.2
    wildtype_beta = wildtype_betas[scenario]
    init_I = initial_weekly_wildtype_incidences[scenario] * popsize

    pl = plotSeqrateVsTravel(variant_beta, wildtype_beta, init_I, 0.05)
    plots[[idx]] = pl[[1]] + 
      ggtitle(bquote(italic(R[e]) ~ " = " ~ .(variant_r0s[variant_r0_idx]))) + 
      theme(plot.title=element_text(size=8,hjust=0.5))
    dfs_travel[[idx]] = pl[[2]]
    dfs_seq[[idx]] = pl[[3]]
    idx = idx + 1
  #}
}

wrap_plots(plots,nrow=2,ncol=2) + plot_layout(guides='collect') & theme(legend.position='bottom')
#ggsave("SuppFig_9.pdf",width=2080,height=1500,units='px',dpi=320)

# Collate all detection times through local sequencing
seq_detection_time = do.call(rbind,dfs_seq)
seq_detection_time[seq_detection_time$variant_beta==0.4,] #For R0 = 2
seq_detection_time[seq_detection_time$variant_beta==0.26,] #For r0 = 1.3


seq_df = do.call(rbind,dfs_travel)
seq_df = seq_df[seq_df$country==T,]

# Get distribution of travel detection times for alpha = 0.1 and R0 = 2
med_2 = median(seq_df[seq_df$psamp==0.1 & seq_df$variant_beta==0.4,]$day,na.rm=T)
iqr_2 = IQR(seq_df[seq_df$psamp==0.1 & seq_df$variant_beta==0.4,]$day,na.rm=T)
c(med_2-iqr_2,med_2,med_2+iqr_2)

# Get distribution of travel detection times for alpha = 0.01 and R0 = 2
med_2 = median(seq_df[seq_df$psamp==0.01 & seq_df$variant_beta==0.4,]$day,na.rm=T)
iqr_2 = IQR(seq_df[seq_df$psamp==0.01 & seq_df$variant_beta==0.4,]$day,na.rm=T)
c(med_2-iqr_2,med_2,med_2+iqr_2)

# Get distribution of travel detection times for alpha = 0.1 and R0 = 1.3
med_2 = median(seq_df[seq_df$psamp==0.1 & seq_df$variant_beta==0.26,]$day,na.rm=T)
iqr_2 = IQR(seq_df[seq_df$psamp==0.1 & seq_df$variant_beta==0.26,]$day,na.rm=T)
c(med_2-iqr_2,med_2,med_2+iqr_2)

# Get distribution of travel detection times for alpha = 0.01 and R0 = 1.3
med_2 = median(seq_df[seq_df$psamp==0.01 & seq_df$variant_beta==0.26,]$day,na.rm=T)
iqr_2 = IQR(seq_df[seq_df$psamp==0.01 & seq_df$variant_beta==0.26,]$day,na.rm=T)
c(med_2-iqr_2,med_2,med_2+iqr_2)


    