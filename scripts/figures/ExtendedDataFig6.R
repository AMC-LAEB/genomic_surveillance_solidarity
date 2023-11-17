library(ggplot2)
library(lemon)
library(scales)
library(ggplot2)
library(dplyr)
library(ggsci)
library(patchwork)


### Read country information
country_data = readRDS("country_data.rds")
cd = country_data[[1]]
continents = cd$Continent
countries = cd$Country

#Plot  theme
width_hi = 1
width_lo = 0.3

fillcol = 'grey95'
colors = c(RColorBrewer::brewer.pal(7,"Dark2")[c(1,2,3,7)],RColorBrewer::brewer.pal(3,"Set1")[2])


thm = theme(axis.text = element_text(color="black",size=5),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            panel.grid.major = element_blank(),
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

# Read outputs
detection_df = readRDS("detection_outputs.rds")
detection_df = detection_df[detection_df$travel_rate=="mean",]
detection_df = detection_df[detection_df$phi == 1,]

detection_df$onsetcontinent = continents[detection_df$onsetcountry]

detection_df$scenario = NA
detection_df[detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.001,]$scenario = 1.2
detection_df[detection_df$wildtype_beta==0.21 & detection_df$base_prev==0.002,]$scenario = 1.3
detection_df[detection_df$wildtype_beta==0.22 & detection_df$base_prev==0.005,]$scenario = 1.6
detection_df[detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.02,]$scenario = 2


### Group time to detection by r0 and surveillance strategy
detection_df_grouped = detection_df %>% 
  group_by(variant_r0,strategy,phi,scenario) %>% 
  dplyr::summarize(Lo = quantile(detection_day,0.025,na.rm=T),
                   Hi=quantile(detection_day,0.975,na.rm=T),
                   Mean = mean(detection_day,na.rm=T),
                   Q25=quantile(detection_day,0.25,na.rm=T),
                   Q75=quantile(detection_day,0.75,na.rm=T),
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))
                   
                   
detection_df_grouped$strategy = factor(detection_df_grouped$strategy,levels=c("min0maxInf","min0max30","min2maxInf","min2max30","double"))


detection_df_grouped$scenario = factor(detection_df_grouped$scenario, 
                                       levels=c(1.2,1.3,1.6,2),
                                       labels=c(expression(paste("wt ",italic(R[e])," = 1, wt prevalence = 0.1%")),
                                                expression(paste("wt ",italic(R[e])," = 1.05, wt prevalence = 0.2%")),
                                                expression(paste("wt ",italic(R[e])," = 1.1, wt prevalence = 0.5%")),
                                                expression(paste("wt ",italic(R[e])," = 1, wt prevalence = 2%"))))


time = ggplot(detection_df_grouped)+ 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=330,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=330,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_lo) + thm +
  theme(legend.position = 'top')+ theme(legend.key.width = unit(.2,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection")+ 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(10,0,10,0)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_hi) +
  scale_color_manual(labels=c(1:5),values=colors) + 
  guides(color=guide_legend(nrow=1),fill='none') + theme(legend.key.size=unit(.7,'lines')) + labs(color="Strategy") +theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0)) + 
  facet_rep_grid(.~scenario,labeller = label_parsed,repeat.tick.labels = T)


### Group incidence by detection by r0 and surveillance strategy
detection_df_grouped = detection_df %>% 
  group_by(variant_r0,strategy,phi,scenario) %>% 
  dplyr::summarize(Lo = quantile(inc,0.025,na.rm=T),
                   Hi=quantile(inc,0.975,na.rm=T),
                   Mean = mean(inc,na.rm=T),
                   Q25=quantile(inc,0.25,na.rm=T),
                   Q75=quantile(inc,0.75,na.rm=T),
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))

detection_df_grouped$strategy = factor(detection_df_grouped$strategy,levels=c("min0maxInf","min0max30","min2maxInf","min2max30","double"))


detection_df_grouped$scenario = factor(detection_df_grouped$scenario, 
                                       levels=c(1.2,1.3,1.6,2),
                                       labels=c(expression(paste("wt ",italic(R[e])," = 1, wt prevalence = 0.1%")),
                                                expression(paste("wt ",italic(R[e])," = 1.05, wt prevalence = 0.2%")),
                                                expression(paste("wt ",italic(R[e])," = 1.1, wt prevalence = 0.5%")),
                                                expression(paste("wt ",italic(R[e])," = 1, wt prevalence = 2%"))))



inc = ggplot(detection_df_grouped) + 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_lo) + thm +
  scale_color_manual(values=colors,labels=c("min0maxInf" = 1,"min0max30" = 2, "min2maxInf" = 3,"min2max30" =4,"double"=5)) + 
  theme(legend.position = 'top') + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Global variant infections\nby day of detection") + theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,0)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_hi) + 
  scale_y_log10(breaks=c(10,100,1000,1e4,1e5,1e6,1e7),labels=trans_format("log10", math_format(10^.x)))+ guides(color=guide_legend(nrow=1))+ theme(legend.key.width = unit(.2,"cm")) + labs(color="Strategy") +
  annotation_logticks(side='l',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + facet_rep_grid(.~scenario,labeller = label_parsed,repeat.tick.labels = T)


time/inc

ggsave("figures/ExtendedDataFig6.pdf",width=180,units="mm",height=120)
