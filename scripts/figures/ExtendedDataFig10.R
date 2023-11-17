library(ggplot2)
library(lemon)
library(scales)
library(ggplot2)
library(dplyr)
library(ggsci)

# Plot theme
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
            legend.spacing.y = unit(.001, 'cm'),
            legend.title=element_text(size=5),
            legend.spacing.x = unit(0.05,"cm"),
            legend.box.spacing = unit(0, "pt"),
            legend.key=element_rect(fill="white"),
            legend.margin=margin(-1.6,0,-1.6,0))

# p = circulating proportion
# rho = prevalence
# N = popsize
# r = seqrate
# s = number

computeError <- function(p, rho, N, r, s, alpha=0.05, c=1e6){
  Za = qnorm(1-(alpha/2))
  q = 1 - p
  e =  Za * sqrt((p*q/r) * (c*rho*s - r)/(N*rho*s - 1))
  return(e)
}



# Compute error for varying sequencing rate, true proportion, confidence level

df = data.frame(popsize=NA, seqrate = NA, error = NA, prop = NA, conf = NA)

seqrates = 10^(seq(-3,3,0.1))

for (popsize in c(5e6,50e6,200e6)){
  for (prop in c(0.1,0.5,0.05)){
    for (conf in c(0.05,0.5)){
      seqrates_new = seqrates[(popsize/1e6)*seqrates >= 1]
      new_df = data.frame(popsize=popsize, seqrate = seqrates_new, error = computeError(prop,0.01,popsize,seqrates_new,1,alpha=conf),prop = prop,conf=conf)
      df = rbind(df,new_df)
    }
  }
}


df = df[!(is.na(df$popsize)),]

df$popsize = factor(df$popsize, levels=c(5e6,5e7,2e8),labels=c(expression(paste("N = 5 million")),expression(paste("N = 50 million")),expression(paste("N = 200 million"))))

ggplot(df,aes(x=seqrate,y=error,group=factor(prop))) + 
  geom_line(data=df[df$conf==0.5,],aes(col=factor(prop),lty="1"),size=lwd) + 
  geom_line(data=df[df$conf==0.05,],aes(col=factor(prop),lty="2"),size=lwd) + 
  coord_cartesian(clip = "off") + thm + 
  scale_color_brewer(palette='Set2') + xlab("Sequencing rate (S/M/wk)") + 
  ylab("Error in estimated variant proportion")  + scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  theme(legend.position = 'top') + scale_linetype_manual(labels=c("50%", "95%"),values=c(1,2)) +
  labs(col="Variant proportion", linetype = "Confidence level") + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) +
  annotation_logticks(side='b',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') + facet_rep_wrap(.~popsize,labeller=label_parsed,repeat.tick.labels = T) +
  theme(legend.position = 'top')


ggsave("figures/ExtendedDataFig10.pdf",width=150,units="mm",height=70)
