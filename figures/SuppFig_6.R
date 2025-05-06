library(ggplot2)
library(lemon)
library(scales)
library(latex2exp)


# Get sequencing rate info
country_data = readRDS("../outputs/country_data.rds") 
cd = country_data[[2]]
high_median = median(cd[cd$Income=="High income",]$Seqrate,na.rm=T)
low_median = median(cd[cd$Income=="Low income",]$Seqrate,na.rm=T)


thm = theme(axis.text = element_text(color="black",size=7),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            panel.grid.major = element_line(linewidth=0.1,color='lightgrey'),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill='white'),
            legend.position='none',
            axis.title= element_text(size=7),
            strip.background=element_rect(colour="white",fill="white"),
            strip.text = element_text(size = 5),
            legend.text = element_text(size =7),
            legend.title=element_text(size=7),
            legend.spacing.x = unit(0.05,"cm"),
            legend.key=element_rect(fill="white"))

lwd = 0.3


time_to_detection_df = readRDS("../outputs/time_to_detection.rds")
time_to_detection_df = time_to_detection_df[!(is.na(time_to_detection_df$variant_r0)),]


time_to_detection_df$scenario = NA
time_to_detection_df[time_to_detection_df$wildtype_r0==1 & time_to_detection_df$initial_wildtype_prevalence==0.001,]$scenario = 1.2
time_to_detection_df[time_to_detection_df$wildtype_r0==1.05 & time_to_detection_df$initial_wildtype_prevalence==0.002,]$scenario = 1.3
time_to_detection_df[time_to_detection_df$wildtype_r0==1.1 & time_to_detection_df$initial_wildtype_prevalence==0.005,]$scenario = 1.6
time_to_detection_df[time_to_detection_df$wildtype_r0==1 & time_to_detection_df$initial_wildtype_prevalence==0.02,]$scenario = 2

time_to_detection_df$scenario = factor(time_to_detection_df$scenario, 
                                       levels=c(1.2,1.3,1.6,2),
                                       labels=c(expression(paste("wt ",italic(R[e])," = 1, wt prevalence = 0.1%")),
                                                expression(paste("wt ",italic(R[e])," = 1.05, wt prevalence = 0.2%")),
                                                expression(paste("wt ",italic(R[e])," = 1.1, wt prevalence = 0.5%")),
                                                expression(paste("wt ",italic(R[e])," = 1, wt prevalence = 2%"))))



incidence_at_detection_plot = ggplot(time_to_detection_df,aes(x=seqrate,y=inc_at_detection,group=variant_r0)) + 
  geom_line(aes(col=variant_r0),size=lwd) + 
  coord_cartesian(clip = "off") + 
  thm + scale_color_brewer(palette='RdYlBu',direction=-1) + geom_rug(data=cd[cd$Seqrate>0.01&cd$Seqrate<1000,],aes(x=Seqrate),inherit.aes=F,size=0.1) + xlab("Sequencing rate (S/M/wk)") + 
  theme(legend.position = 'top') + 
  ylab("Variant infections by\nday of detection")  + scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,1e4,1e5,1e6,1e7,1e8),labels=trans_format("log10", math_format(10^.x))) +
  labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) + 
  annotation_logticks(side='lb',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') + #+ theme(legend.box.margin=margin(-5,0,-5,0))
  facet_rep_grid(cols=vars(scenario),labeller = label_parsed,repeat.tick.labels = T) +
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10))



incidence_at_detection_plot_dn = ggplot(time_to_detection_df,aes(x=seqrate,y=inc_at_detection,group=variant_r0)) + 
  geom_line(aes(x=seqrate,y=inc_at_detection_dn,group=variant_r0,col=variant_r0),size=lwd) + 
  coord_cartesian(clip = "off") + guides(color='none') +
  thm + scale_color_brewer(palette='RdYlBu',direction=-1) + geom_rug(data=cd[cd$Seqrate>0.01&cd$Seqrate<1000,],aes(x=Seqrate),inherit.aes=F,size=0.1) + xlab("Sequencing rate (S/M/wk)") + 
  ylab("Reduction in variant infections\nby day of detection per S/M/wk added")  + scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000,1e4,1e5,1e6,1e7,1e8),labels=trans_format("log10", math_format(10^.x))) +
  theme(legend.position = 'none') + scale_linetype_manual(labels=c(TeX("$I(\\tau)$",italic=T),TeX("$\\Delta I(\\tau) / \\Delta n$",italic=TRUE)),values=c(1,2)) +
  labs(col=expression("Variant"~ italic("R"[e]))) + theme(legend.box='vertical') + theme(legend.spacing.y=unit(-.1,'cm')) + 
  annotation_logticks(side='lb',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') + 
  facet_rep_grid(cols=vars(scenario),labeller = label_parsed,repeat.tick.labels = T)


incidence_at_detection_plot+incidence_at_detection_plot_dn+plot_layout(design='A\nB') & plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=8,face='bold'))

ggsave("SuppFig_6.pdf",width=2080,height=1500,units='px',dpi=320)


