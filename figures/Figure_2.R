library(ggplot2)
library(scales)
library(latex2exp)
library(patchwork)

# Get sequencing rate info
country_data = readRDS("../outputs/country_data.rds") 
cd = country_data[[2]]
high_median = median(cd[cd$Income=="High income",]$Seqrate,na.rm=T)
low_median = median(cd[cd$Income=="Low income",]$Seqrate,na.rm=T)

# Read single-country simulation output
time_to_detection_df = readRDS("../outputs/time_to_detection.RDS")
effective_increase_df = readRDS("../outputs/effective_increase.RDS")

time_to_detection_df = time_to_detection_df[!(is.na(time_to_detection_df$variant_r0)),]
effective_increase_df = effective_increase_df[!(is.na(effective_increase_df$variant_r0)),]

# Keep only scenarios described in main text

time_to_detection_df$scenario = NA
time_to_detection_df[time_to_detection_df$variant_r0==1.2 & time_to_detection_df$wildtype_r0==1 & time_to_detection_df$initial_wildtype_prevalence==0.001,]$scenario = 1.2
time_to_detection_df[time_to_detection_df$variant_r0==1.3 & time_to_detection_df$wildtype_r0==1.05 & time_to_detection_df$initial_wildtype_prevalence==0.002,]$scenario = 1.3
time_to_detection_df[time_to_detection_df$variant_r0==1.6 & time_to_detection_df$wildtype_r0==1.1 & time_to_detection_df$initial_wildtype_prevalence==0.005,]$scenario = 1.6
time_to_detection_df[time_to_detection_df$variant_r0==2 & time_to_detection_df$wildtype_r0==1 & time_to_detection_df$initial_wildtype_prevalence==0.02,]$scenario = 2

effective_increase_df$scenario = NA
effective_increase_df[effective_increase_df$variant_r0==1.2 & effective_increase_df$wildtype_r0==1 & effective_increase_df$initial_wildtype_prevalence==0.001,]$scenario = 1.2
effective_increase_df[effective_increase_df$variant_r0==1.3 & effective_increase_df$wildtype_r0==1.05 & effective_increase_df$initial_wildtype_prevalence==0.002,]$scenario = 1.3
effective_increase_df[effective_increase_df$variant_r0==1.6 & effective_increase_df$wildtype_r0==1.1 & effective_increase_df$initial_wildtype_prevalence==0.005,]$scenario = 1.6
effective_increase_df[effective_increase_df$variant_r0==2 & effective_increase_df$wildtype_r0==1 & effective_increase_df$initial_wildtype_prevalence==0.02,]$scenario = 2

time_to_detection_df = time_to_detection_df[!(is.na(time_to_detection_df$scenario)),]
effective_increase_df = effective_increase_df[!(is.na(effective_increase_df$scenario)),]


### Plot theme

thm =     theme(axis.text = element_text(color="black",size=7),
                panel.spacing.x = unit(1, "mm"),
                axis.line=element_line(color='black',size = 0.2),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                legend.position='none',
                axis.title= element_text(size=7),
                strip.background=element_rect(colour="black",fill="white"),
                strip.text.x = element_text(size = 7),
                legend.text=element_text(size=7),
                legend.title =element_text(size=7),
                legend.key=element_rect(fill="white"),
                legend.key.width = unit(.5, "line"),
                legend.spacing.y = unit(.001, 'cm'),
                legend.box.spacing = unit(0, "pt"),
                legend.margin=margin(-1.6,0,-1.6,0))

lwd = 0.5
lwd_line = 0.3

### Fig. 2a
time_to_detection_plot = ggplot(time_to_detection_df,aes(x=seqrate,y=day_of_detection,group=variant_r0)) + 
  geom_line(aes(col=variant_r0),size=lwd) +  thm + #coord_cartesian(clip = "off") +
  scale_color_brewer(palette='RdYlBu',direction=-1) + geom_rug(data=cd[cd$Seqrate>0.01 & cd$Seqrate<1000,],aes(x=Seqrate),inherit.aes=F,size=0.1) + xlab("Sequencing rate (S/M/wk)") + 
  ylab("Time to detection (d)")  + scale_x_continuous(breaks=c(0.01,2,10,25,30,50,75,100),limits=c(0.01,100),labels=c(0.01,'\n2',10,25,'\n30',50,75,100)) + 
  theme(legend.position = 'top') +
  labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) +
  guides(color='none') +
  geom_vline(aes(xintercept=low_median),lty=2,linewidth=lwd) + geom_vline(aes(xintercept=high_median),lty=2,linewidth=lwd) + 
  scale_y_continuous(breaks=seq(0,1000,50)) + 
  theme(axis.text.x = element_text(color=c('black','grey50','black','grey50','black','black','black')),
        axis.ticks.x = element_line(color=c('black','grey50','black','grey50','black','black','black'),
                                    size=unit(c(.5,0.2,.5,.2,.5,.5,.5),'cm'))) + 
  geom_vline(xintercept = c(2,30),linewidth=0.1,col='grey50',lty=1) 
y_coord = layer_scales(time_to_detection_plot)$y$range$range[2]
time_to_detection_plot = time_to_detection_plot + annotate("text",x=low_median,y=y_coord*0.95,label="LIC",angle=90,size=7/.pt,vjust=2) +
  annotate("text",x=high_median,y=y_coord*0.95,label="HIC",angle=90,size=7/.pt,vjust=2) 
  



time_to_detection_plot = ggplot(time_to_detection_df,aes(x=seqrate,y=day_of_detection,group=variant_r0)) + 
  geom_line(aes(col=variant_r0),size=lwd) +  thm + #coord_cartesian(clip = "off") +
  scale_color_brewer(palette='RdYlBu',direction=-1) + geom_rug(data=cd[cd$Seqrate>0.01 & cd$Seqrate<1000,],aes(x=Seqrate),inherit.aes=F,size=0.1) + xlab("Sequencing rate (S/M/wk)") + 
  ylab("Time to detection (d)")  + scale_x_continuous(breaks=c(2,10,30,50,75,100),limits=c(0.01,100),labels=c(2,10,30,50,75,100)) + 
  theme(legend.position = 'top') +
  labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) +
  guides(color='none') +
  geom_vline(aes(xintercept=low_median),lty=2,linewidth=lwd_line) + geom_vline(aes(xintercept=high_median),lty=2,linewidth=lwd_line) + 
  scale_y_continuous(breaks=seq(0,1000,50)) + 
  theme(axis.text.x = element_text(color='black'),
        axis.ticks.x = element_line(color='black',
                                    size=unit(c(.5,0.2,.5,.2,.5,.5,.5),'cm'))) + 
  geom_vline(xintercept = c(2,30),linewidth=0.1,col='grey50',lty=1) 
y_coord = layer_scales(time_to_detection_plot)$y$range$range
y_coord = y_coord[1] + 0.95 * diff(y_coord)
time_to_detection_plot = time_to_detection_plot + annotate("text",x=low_median,y=y_coord,label="LIC",angle=90,size=7/.pt,vjust=2) +
  annotate("text",x=high_median,y=y_coord,label="HIC",angle=90,size=7/.pt,vjust=2) 




time_to_detection_dn_plot = ggplot(time_to_detection_df,aes(x=seqrate,y=day_of_detection,group=variant_r0)) + 
  geom_line(aes(x=seqrate,y=day_of_detection_dn,group=variant_r0,col=variant_r0),size=lwd) + 
  coord_cartesian(clip = "off")+
  thm + scale_color_brewer(palette='RdYlBu',direction=-1) + geom_rug(data=cd[cd$Seqrate>0.01&cd$Seqrate<1000,],aes(x=Seqrate),inherit.aes=F,size=0.1) + xlab("Sequencing rate (S/M/wk)") + 
  ylab("Reduction in time to detection\nper S/M/wk added (d)")  + scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  theme(legend.position = 'none') + 
  labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) +
  guides(color='none') + 
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100), labels=trans_format("log10", math_format(10^.x))) + 
  geom_vline(aes(xintercept=low_median),lty=2,linewidth=lwd_line) + geom_vline(aes(xintercept=high_median),lty=2,linewidth=lwd_line) + 
  annotation_logticks(side='bl',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') 
y_coord = layer_scales(time_to_detection_dn_plot)$y$range$range
y_coord = y_coord[1] + 0.95 * diff(y_coord)
time_to_detection_dn_plot = time_to_detection_dn_plot + annotate("text",x=low_median,y=10^(y_coord),label="LIC",angle=90,size=7/.pt,vjust=2) +
  annotate("text",x=high_median,y=10^(y_coord),label="HIC",angle=90,size=7/.pt,vjust=2) 



### Fig. 2b
incidence_at_detection_plot = ggplot(time_to_detection_df,aes(x=seqrate,y=inc_at_detection,group=variant_r0)) + 
  geom_line(aes(col=variant_r0),size=lwd) + 
  coord_cartesian(clip = "off") + guides(color='none') + 
  thm + scale_color_brewer(palette='RdYlBu',direction=-1) + geom_rug(data=cd[cd$Seqrate>0.01&cd$Seqrate<1000,],aes(x=Seqrate),inherit.aes=F,size=0.1) + xlab("Sequencing rate (S/M/wk)") + 
  ylab("Variant infections by detection day")  + scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,1e4,1e5,1e6,1e7,1e8),labels=trans_format("log10", math_format(10^.x))) +
  theme(legend.position = 'none') + 
  labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) + 
  geom_vline(aes(xintercept=low_median),lty=2,linewidth=lwd_line) + geom_vline(aes(xintercept=high_median),lty=2,linewidth=lwd_line) + 
  annotation_logticks(side='lb',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') 
y_coord = layer_scales(incidence_at_detection_plot)$y$range$range
y_coord = y_coord[1] + 0.95 * diff(y_coord)
incidence_at_detection_plot = incidence_at_detection_plot + annotate("text",x=low_median,y=10^(y_coord),label="LIC",angle=90,size=7/.pt,vjust=2) +
 annotate("text",x=high_median,y=10^(y_coord),label="HIC",angle=90,size=7/.pt,vjust=2) 
  

incidence_at_detection_plot_dn = ggplot(time_to_detection_df,aes(x=seqrate,y=inc_at_detection,group=variant_r0)) + 
  geom_line(aes(x=seqrate,y=inc_at_detection_dn,group=variant_r0,col=variant_r0),size=lwd) + 
  coord_cartesian(clip = "off") + theme(legend.position='none') +#guides(color='none') +
  thm + scale_color_brewer(palette='RdYlBu',direction=-1) + geom_rug(data=cd[cd$Seqrate>0.01&cd$Seqrate<1000,],aes(x=Seqrate),inherit.aes=F,size=0.1) + xlab("Sequencing rate (S/M/wk)") + 
  ylab("Reduction in variant infections by day\nof detection per S/M/wk added")  + scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,1e4,1e5,1e6,1e7,1e8),labels=trans_format("log10", math_format(10^.x))) +
  theme(legend.position = 'none') + 
  labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) + 
  geom_vline(aes(xintercept=low_median),lty=2,linewidth=lwd_line) + geom_vline(aes(xintercept=high_median),lty=2,linewidth=lwd_line) + 
  annotation_logticks(side='lb',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') + guides(color=guide_legend(override.aes = list(linewidth=.8)))
y_coord = layer_scales(incidence_at_detection_plot_dn)$y$range$range
y_coord = y_coord[1] + 0.95 * diff(y_coord)
incidence_at_detection_plot_dn = incidence_at_detection_plot_dn + annotate("text",x=low_median,y=10^(y_coord),label="LIC",angle=90,size=7/.pt,vjust=2) + 
annotate("text",x=high_median,y=10^(y_coord),label="HIC",angle=90,size=7/.pt,vjust=2) 


### Fig. 2c
effective_increase_plot = ggplot(effective_increase_df,aes(x=tat_reduction,y=equivalent_increase,group=factor(variant_r0))) + 
  geom_line(aes(color=factor(variant_r0)),size=lwd) + thm + theme(legend.position = 'none') + scale_color_brewer(palette='RdYlBu',direction=-1) + 
  xlab("Reduction in turnaround time (d)") + ylab("Equivalent fold increase\nin sequencing rate") +  
  scale_y_log10(breaks=c(1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) + labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.position = 'none') + 
  annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') 

# time_to_detection_plot+time_to_detection_dn_plot+incidence_at_detection_plot+incidence_at_detection_plot_dn + effective_increase_plot +
#   plot_layout(guides='collect',nrow=1) & theme(legend.position = 'bottom',legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10)) &
#   plot_annotation(theme = theme(plot.margin = unit(c(0,0,0,0),'cm')))

### Compute the equivalent increase in sequencing rate
print("Equivalent fold increase in sequencing rate")
print(effective_increase_df[effective_increase_df$scenario==2 & effective_increase_df$tat_reduction==21,]$equivalent_increase)
print(effective_increase_df[effective_increase_df$scenario==1.2 & effective_increase_df$tat_reduction==21,]$equivalent_increase)



time_to_detection_plot+theme(plot.margin=unit(c(0,0.1,0,-0.1),"cm"))+
  time_to_detection_dn_plot+theme(plot.margin=unit(c(0,0.1,0,-0.1),"cm"))+
  incidence_at_detection_plot+theme(plot.margin=unit(c(0,0.1,0,-0.1),"cm"))+
  effective_increase_plot+theme(plot.margin=unit(c(0,0.15,0,-0.1),"cm"))+
  plot_layout(guides='collect',nrow=1,ncol=4) & theme(legend.position = 'bottom',
                                                      legend.margin=margin(0,0,0,0),
                                                      legend.box.margin=margin(-10,-10,-10,-10)) #&
  #plot_annotation(theme = theme(plot.margin = unit(c(0,-0.1,0,-0.1),'cm')))


ggsave("Figure_2.pdf",width=2080,height=650,units='px',dpi=320)




### Compute the reduction in time to detection for R_e = 1.6 in a LIC and HIC
lic_diff = time_to_detection_df[time_to_detection_df$seqrate==low_median & time_to_detection_df$variant_r0==1.6,]$day_of_detection_dn
hic_diff = time_to_detection_df[time_to_detection_df$seqrate==high_median & time_to_detection_df$variant_r0==1.6,]$day_of_detection_dn

print("Reduction in time to detection for high-income country")
print(hic_diff)
print("Reduction in time to detection for low-income country")
print(lic_diff)

### Compute the incidence by time of detection for R_e = 1.6 in a LIC and HIC
lic_inc = time_to_detection_df[time_to_detection_df$seqrate==low_median & time_to_detection_df$variant_r0==1.6,]$inc_at_detection
hic_inc = time_to_detection_df[time_to_detection_df$seqrate==high_median & time_to_detection_df$variant_r0==1.6,]$inc_at_detection

print("Incidence by detection for high-income country")
print(hic_inc)
print("Incidence by deection for low-income country")
print(lic_inc)

### Compute the reduction in time to detection for R_e = 1.6 in a LIC and HIC
lic_inc_diff = time_to_detection_df[time_to_detection_df$seqrate==low_median & time_to_detection_df$variant_r0==1.6,]$inc_at_detection_dn
hic_inc_diff = time_to_detection_df[time_to_detection_df$seqrate==high_median & time_to_detection_df$variant_r0==1.6,]$inc_at_detection_dn

print("Reduction in incidence by detection for high-income country")
print(hic_inc_diff)
print(lic_inc_diff)

