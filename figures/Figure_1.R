library(ggridges)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggsci)

#Get country data

country_data = readRDS("../outputs/country_data.rds")
cd = country_data[[1]]
continents = cd$Continent
countries = cd$Country

#Plotting params

width_hi = 1
width_lo = 0.3

fillcol = 'grey95'
nday = 1000

f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
cols = f('Set1')[c(1:5,7)]


thm = theme(axis.text = element_text(color="black",size=6),
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
            legend.text = element_text(size = 6),
            legend.title=element_blank(),
            legend.spacing.x = unit(0.05,"cm"),
            legend.key=element_rect(fill="white"))


### Read genomic surveillance outputs
detection_df = readRDS("../outputs/combined_outputs.rds")
detection_df[detection_df$detection_day==0,]$detection_day = NA
detection_df = detection_df[detection_df$travel_rate=="mean" & detection_df$phi==1,]

# Keep only scenarios described in the main text
detection_df$scenario = NA
detection_df[detection_df$variant_r0==1.2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.001,]$scenario = 1.2
detection_df[detection_df$variant_r0==1.3 & detection_df$wildtype_beta==0.21 & detection_df$base_prev==0.002,]$scenario = 1.3
detection_df[detection_df$variant_r0==1.6 & detection_df$wildtype_beta==0.22 & detection_df$base_prev==0.005,]$scenario = 1.6
detection_df[detection_df$variant_r0==2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.02,]$scenario = 2

detection_df = detection_df[!(is.na(detection_df$scenario)),]
detection_df = detection_df[!(is.na(detection_df$detection_day)),]

detection_df$onsetcontinent = continents[detection_df$onset_country]



### Plot time to detection ridgeplots
time_to_detection = ggplot(detection_df[detection_df$strategy=='min0maxInf' & detection_df$phi == 1,],aes(y=factor(variant_r0),fill=stat(x))) + 
  geom_density_ridges_gradient(rel_min_height=0.0000001,quantile_lines=TRUE,quantiles=c(0.025,0.5,0.975),aes(x=detection_day,group=factor(variant_r0)), 
                               vline_size=0.2, alpha = .3, color = 'black', size=0.1,fill='grey90') +
  thm + xlab("Day of detection") + ylab(expression("Variant"~ italic("R"[e]))) + xlim(c(0,300)) + theme(plot.margin = unit(c(.2,.1,.2,.1),'cm'))

### Plot incidence by detection ridgeplots
incidence_at_detection = ggplot(detection_df[detection_df$strategy=='min0maxInf' & detection_df$phi == 1,],aes(y=factor(variant_r0),fill=stat(x))) + 
  geom_density_ridges_gradient(rel_min_height=0.001,quantile_lines=TRUE,quantiles=c(0.025,0.5,0.975),aes(x=detection_infections,group=factor(variant_r0)),
                               vline_size=0.2, alpha = .3, color = 'black', size=0.1,fill='grey90') +thm +
  xlab("Global variant infections\nby day of detection") + ylab(expression("Variant"~ italic("R"[e]))) + scale_x_log10(breaks=10^seq(0,8,2),labels=trans_format("log10", math_format(10^.x))) + 
  theme(plot.margin = unit(c(.2,.1,.2,.1),'cm')) + annotation_logticks(side='b',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) +
  coord_cartesian(clip = "off")


### Group outputs by continent for time to detection and compute probability that it is detected in its origin continent
detection_df_grouped = detection_df %>% 
  group_by(onsetcontinent,variant_r0,strategy) %>% 
  dplyr::summarize(Lo = quantile(detection_day,0.025,na.rm=T),
                   Hi=quantile(detection_day,0.975,na.rm=T),
                   Q25=quantile(detection_day,0.25,na.rm=T),
                   Q75=quantile(detection_day,0.75,na.rm=T),
                   Mean = mean(detection_day,na.rm=T),
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))

### Plot time to detection by continent 
time_to_detection_per_continent = ggplot(detection_df_grouped[detection_df_grouped$strategy=="min0maxInf",])+ 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=1),linewidth=width_lo) + thm +
  scale_color_manual(values=cols,labels=c("AF","AS","EU","NA","OC","SA")) + scale_fill_manual(values=cols,guide='none') + 
  theme(legend.position = 'top') + theme(legend.key.width = unit(.15,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection") + guides(color=guide_legend(nrow=1,override.aes = list(shape = NA)),fill='none') + 
  theme(legend.margin=margin(-10,0,-10,0),legend.box.margin=margin(-10,-20,-10,-20)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=1),linewidth=width_hi) +
  theme(plot.margin = unit(c(0,.1,0,.1),'cm')) +
  geom_point(aes(x=factor(variant_r0),y=Mean,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=1),cex=1,pch=21,fill='white') 


### Plot probability of local detection per continent
p_local_detection_per_continent = ggplot(detection_df_grouped[detection_df_grouped$strategy=="min0maxInf",]) + 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) + thm +
  geom_point(aes(x=factor(variant_r0),y=LocalDetect_Continent,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=.7),pch=3,size=1) +
  scale_color_manual(values=cols,labels=c("AF","AS","EU","NA","OC","SA")) +
  theme(legend.position = 'top')+ theme(legend.key.width = unit(.15,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Probability of first detection\nin origin continent")+ guides(color=guide_legend(nrow=1)) + 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-20,-10,-40)) + ylim(c(0,1))

### Group outputs by continent for incidence at detection
detection_df_grouped = detection_df %>% 
  group_by(onsetcontinent,variant_r0,strategy) %>% 
  dplyr::summarize(Lo = quantile(detection_infections,0.025,na.rm=T),
                   Mean=mean(detection_infections,na.rm=T),
                   Hi=quantile(detection_infections,0.975,na.rm=T),
                   Q25=quantile(detection_infections,0.25,na.rm=T),
                   Q75=quantile(detection_infections,0.75,na.rm=T))

### Plot incidence by detection per continent
incidence_at_detection_per_continent = ggplot(detection_df_grouped[detection_df_grouped$strategy=="min0maxInf",])+ 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=1),linewidth=width_lo) + thm +
  scale_color_manual(values=cols,labels=c("AF","AS","EU","NA","OC","SA")) + scale_fill_manual(values=cols,guide='none') + 
  theme(legend.position = 'top') + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Global variant infections\nby day of detection") + theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,-20)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=1),linewidth=width_hi) + 
  scale_y_log10(breaks=c(1,10,100,1000,1e4,1e5,1e6,1e7),labels=trans_format("log10", math_format(10^.x)))+ guides(color=guide_legend(nrow=1,override.aes = list(shape = NA)),fill='none')+ theme(legend.key.width = unit(.15,"cm")) + 
  annotation_logticks(side='l',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) +
  geom_point(aes(x=factor(variant_r0),y=Mean,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=1),cex=1,pch=21,fill='white') 



### Group outputs by onset country
detection_df_grouped_country = detection_df[detection_df$strategy=="min0maxInf",] %>% 
  group_by(onset_country,variant_r0,onsetcontinent) %>% 
  dplyr::summarize(Mean = mean(detection_day,na.rm=T),
                   Inc = mean(detection_infections,na.rm=T),
                   Low=quantile(detection_day,0.95,na.rm=T),
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))

detection_df_grouped_country$seqrate = cd$Seqrate[detection_df_grouped_country$onset_country]
detection_df_grouped_country$seqrate[detection_df_grouped_country$seqrate==0] = NA

### Plot time to detection by origin country's sequencing rate
time_to_detection_by_seqrate = ggplot(detection_df_grouped_country) + 
  labs(color=expression("Variant"~ italic("R"[e]))) + 
  geom_point(aes(x=seqrate,y=Mean,fill=factor(variant_r0)),pch=21,stroke=0.1,alpha=0.8,col='black',cex=.7) + theme(legend.position = 'none')  +
  thm + theme(legend.position='top') + xlab("Origin country sequencing\nrate (S/M/wk)") + ylab("Mean day of detection") +
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) + coord_cartesian(clip = "off") +  scale_size(range=c(0.01,4)) + 
  geom_smooth(aes(x=seqrate,y=Mean,group=factor(variant_r0),color=factor(variant_r0)),method='loess',se=F,linewidth=.5,lty=2,alpha=0.3) + scale_color_brewer(palette='RdYlBu',direction=-1) + 
  scale_fill_brewer(palette='RdYlBu',direction=-1) + 
  guides(fill='none',col = guide_legend(barheight=.7)) +
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-20,-10,-30),legend.title = element_text(size=5),legend.key.width = unit(.2,"cm")) + ylim(c(0,250)) + 
  annotation_logticks(side='b',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1)+ guides(color=guide_legend(override.aes = list(linewidth=.8)))



cd = country_data[[2]]
cd = cd[!(is.na(cd$Population)),]

#Plot theme
thm1 =     theme(axis.text = element_text(color="black",size=6),
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
                 legend.text=element_text(size=6),
                 legend.title =element_text(size=7),
                 legend.key=element_rect(fill="white"),
                 legend.key.width = unit(.5, "line"),
                 legend.spacing.y = unit(.001, 'cm'),
                 legend.box.spacing = unit(0, "pt"),
                 legend.margin=margin(1,1,1,1))

lty = 1
lw = 0.2
lcol='lightgrey'

# Plot distribution of sequencing rates
cd$Seqrate_2 = cd$Seqrate
cd[cd$Seqrate_2==0,]$Seqrate_2=1e-4
sequencing_rate_histogram = ggplot(cd,aes(x=Seqrate_2,group=Continent,fill=Continent)) + geom_histogram(linewidth=0.1,col='black',bins=15,alpha=0.7) + thm1  +
  scale_x_log10(breaks=c(0.0001,0.01,0.1,1,10,100,1000),labels=c(0,expression(10^-2),expression(10^-1),expression(10^-0),expression(10^1),expression(10^2),expression(10^3))) + scale_fill_manual(values=cols,labels=c("AF","AS","EU","NA","OC","SA")) + xlab("Sequencing rate (S/M/wk)") + ylab("Count")+
  theme(legend.position='top') + theme(legend.title=element_blank()) + theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm")) +
  coord_cartesian(clip='off') + geom_vline(xintercept=5e-4,linewidth=0.2,lty=2)

# Plot distribution of turnaround times
turnaround_time_histogram = ggplot(cd[!(is.na(cd$Median_TAT)),],aes(x=Median_TAT,group=Continent,fill=Continent)) + geom_histogram(linewidth=0.1,col='black',bins=15,alpha=0.7) + thm1 +
  scale_fill_manual(values=cols,labels=c("AF","AS","EU","NA","OC","SA")) + 
  xlab("Median turnaround time (d)") + ylab("Count")  + theme(legend.position='top') + theme(legend.title=element_blank()) +
  theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm")) + scale_x_log10(breaks=c(1,10,30,100,300),labels=c(1,10,30,100,300)) + 
  annotation_logticks(side='b',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1)+ coord_cartesian(clip='off')


ggarrange(ggarrange(sequencing_rate_histogram,turnaround_time_histogram,((time_to_detection|incidence_at_detection) + plot_annotation(theme = theme(plot.margin = unit(c(0,.1,-.2,.1),'cm')))),
                    widths=c(0.22,0.22,0.4,0.4),nrow=1),
          ggarrange((p_local_detection_per_continent|time_to_detection_per_continent|incidence_at_detection_per_continent|time_to_detection_by_seqrate)+plot_layout(widths=c(0.3,0.4,0.4,0.3)),nrow=1),ncol=1,heights=c(0.45,0.55))


ggsave("Figure_1.pdf",width=2080,height=1200,units='px',dpi=320)




