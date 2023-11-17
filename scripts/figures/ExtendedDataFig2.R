library(ggplot2)

### Get input parameters (variant Re, wildtype beta, wildtype prevalence)

v_1.2 = c(1.2, 0.2, 0.001)
v_1.3 = c(1.3, 0.21, 0.002)
v_1.6 = c(1.6, 0.22, 0.005)
v_2 = c(2, 0.2, 0.02)

popsize = 100e6

params = list(v_1.2, v_1.3, v_1.6, v_2)

### Simulates the SIR epidemic dynamics of two viruses

simulateEpidemic <- function(state_vec, params, ndays,N, introductionday){
  
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


### Simulate for all four scenarios in main text
df = data.frame(variant_r0=NA,t=NA,wt=NA,mt=NA)

for (i in 1:4){
  init_I = params[[i]][3] * popsize
  wildtype_beta = params[[i]][2]
  variant_beta = params[[i]][1]*0.2
  
  state_vec = c(popsize-init_I,init_I,0,popsize,1,0)
  epi_params = c(wildtype_beta, variant_beta, 0.2)
  
  sim = simulateEpidemic(state_vec,epi_params,1000,popsize)
  df = rbind(df,data.frame(variant_r0=variant_beta,t=1:1000,wt=sim[[2]][,1],mt=sim[[2]][,2]))
}

df = gather(df,virus,inc,wt:mt)
df$inc = (df$inc/popsize) * 100



### Plot epidemic dynamics
df$variant_r0 = factor(df$variant_r0, levels=c(1.2,1.3,1.6,2)*0.2,labels=c(expression(atop(paste("Variant ",italic(R[e])," = 1.2"), paste("(Wildtype ",italic(R[e])," = 1, wildtype prevalence = 0.1%)"))),
                                                                           expression(atop(paste("Variant ",italic(R[e])," = 1.3"), paste("(Wildtype ",italic(R[e])," = 1.05, wildtype prevalence = 0.2%)"))),
                                                                           expression(atop(paste("Variant ",italic(R[e])," = 1.6"), paste("(Wildtype ",italic(R[e])," = 1.1, wildtype prevalence = 0.5%)"))),
                                                                           expression(atop(paste("Variant ",italic(R[e])," = 2"), paste("(Wildtype ",italic(R[e])," = 1, wildtype prevalence = 2%)")))))

df = df[!(is.na(df$variant_r0)),]
df$virus = factor(df$virus,levels=c("mt","wt"),labels=c("Variant","Wildtype"))

thm = theme(axis.text = element_text(color="black",size=5),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            panel.grid.major = element_line(linewidth=0.1,color='lightgrey'),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill='white'),
            legend.position='none',
            axis.title= element_text(size=5),
            strip.background=element_rect(colour="white",fill="white"),
            strip.text = element_text(size = 5),
            legend.text = element_text(size = 5),
            legend.title=element_text(size=5),
            legend.spacing.x = unit(0.05,"cm"),
            legend.key=element_rect(fill="white"))



ggplot(df,aes(x=t,y=inc,group=virus,col=virus)) + geom_line(linewidth=0.3) + facet_wrap(.~variant_r0,scales='free',labeller=label_parsed) + xlim(c(0,300)) + thm + ylab("Incidence (%)") + 
  xlab("Day") + theme(legend.position='top',legend.title=element_blank()) + theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10)) + scale_color_brewer(palette="Set2")

ggsave("figures/ExtendedDataFig2.pdf",width=100,units="mm",height=80)






