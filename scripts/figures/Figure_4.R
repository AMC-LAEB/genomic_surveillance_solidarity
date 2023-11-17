library(ggplot2)
library(ggsci)
library(latex2exp)
library(lemon)
library(dplyr)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)

# Get country data

country_data = readRDS("country_data.rds")
cd = country_data[[1]]
continents = cd$Continent
countries = cd$Country



# Plotting parameters

thm = theme(axis.text = element_text(color="black",size=5),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position='none',
            axis.title= element_text(size=5),
            strip.background=element_rect(colour="black",fill="white"),
            strip.text.x = element_text(size = 5),
            legend.text = element_text(size = 5),
            legend.title=element_text(size=5),
            legend.spacing.x = unit(0.05,"cm"),
            legend.key=element_rect(fill="white"))

colors = c(RColorBrewer::brewer.pal(7,"Dark2")[c(1,2,3,7)],RColorBrewer::brewer.pal(3,"Set1")[2])
fillcol = 'grey95'
nday = 1000

width_hi = 1
width_lo = 0.3

### Read genomic surveillance outputs
detection_df = readRDS("detection_outputs.rds")
detection_df = detection_df[detection_df$travel_rate=="mean",]
detection_df = detection_df[detection_df$phi == 1,]

detection_df$strategy = factor(detection_df$strategy,levels=c("min0maxInf","min0max30","min2maxInf","min2max30","double"))


# Keep only scenarios described in the main text
detection_df$scenario = NA
detection_df[detection_df$variant_r0==1.2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.001,]$scenario = 1.2
detection_df[detection_df$variant_r0==1.3 & detection_df$wildtype_beta==0.21 & detection_df$base_prev==0.002,]$scenario = 1.3
detection_df[detection_df$variant_r0==1.6 & detection_df$wildtype_beta==0.22 & detection_df$base_prev==0.005,]$scenario = 1.6
detection_df[detection_df$variant_r0==2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.02,]$scenario = 2
detection_df = detection_df[!(is.na(detection_df$scenario)),]
detection_df = detection_df[!(is.na(detection_df$detection_day)),]

detection_df$onsetcontinent = continents[detection_df$onsetcountry]


### Group time to detection by r0 and surveillance strategy
detection_df_grouped = detection_df %>% group_by(variant_r0,strategy) %>% dplyr::summarize(Lo = quantile(detection_day,0.025,na.rm=T),
                                                                                           Hi=quantile(detection_day,0.975,na.rm=T),
                                                                                           Mean = mean(detection_day,na.rm=T),
                                                                                           Q25=quantile(detection_day,0.25,na.rm=T),
                                                                                           Q75=quantile(detection_day,0.75,na.rm=T),
                                                                                           LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))

detection_df_grouped$strategy = factor(detection_df_grouped$strategy,levels=c("min0maxInf","min0max30","min2maxInf","min2max30","double"))


time_to_detection_by_strategy = ggplot(detection_df_grouped) + 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=330,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=330,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_lo) + thm +
  scale_color_manual(values=colors, labels=1:5) + 
  theme(legend.position = 'top') + theme(legend.key.width = unit(.2,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection") + 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(10,0,10,0)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_hi) +
  guides(color=guide_legend(nrow=1),fill='none') + theme(legend.key.size=unit(width_hi,'lines')) + labs(color="Strategy") +theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,0))

### Get legend
legplot = ggplot(detection_df_grouped)+ 
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_hi) + thm +
  scale_color_manual(values=colors,labels=c("1: 2022 baseline","2: 2022 baseline, capped at 30 S/M/wk","3: 2022 baseline + minimum global capacity",  "4: 2022 baseline, capped at 30 S/M/wk + minimum global capacity","5: 2022 baseline, doubled")) +
  guides(color=guide_legend(nrow=1),fill='none') + theme(legend.key.size=unit(1,'lines')) + labs(color="Strategy") + theme(legend.position = 'top')+ theme(legend.key.width = unit(.2,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection")+ 
  theme(legend.margin=margin(-10,0,0,0),legend.box.margin=margin(-10,0,-10,0))

### Plot probability of local detection by strategy by r0
p_local_detection_by_strategy = ggplot(detection_df_grouped)+ 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_point(aes(x=factor(variant_r0),y=LocalDetect_Continent,group=strategy,col=strategy),position=position_dodge(width=.7),pch=3,size=1) + thm +
  scale_color_manual(values=colors,labels=c("min0maxInf" = 1,"min0max30" = 2,"min2maxInf" = 3, "min2max30" =4,"double"=5)) + 
  theme(legend.position = 'top')+ theme(legend.key.width = unit(.2,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Probability of first detection\nin origin continent")+ 
  guides(color=guide_legend(nrow=1)) + 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,10,-10,0)) + labs(color="Strategy") + ylim(c(0.6,1)) 


### Group incidence by r0 and strategy
detection_df_grouped = detection_df %>% group_by(variant_r0,strategy) %>% dplyr::summarize(Lo = quantile(inc,0.025,na.rm=T),
                                                                                           Hi=quantile(inc,0.975,na.rm=T),
                                                                                           Med = quantile(inc,0.5,na.rm=T),
                                                                                           Q25=quantile(inc,0.25,na.rm=T),
                                                                                           Q75=quantile(inc,0.75,na.rm=T),
                                                                                           Mean=mean(inc,na.rm=T))

# Plot incidence at detection by strategy
detection_df_grouped$strategy = factor(detection_df_grouped$strategy,levels=c("min0maxInf","min0max30","min2maxInf","min2max30","double"))
incidence_at_detection_by_strategy = ggplot(detection_df_grouped) + 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_lo) + thm +
  scale_color_manual(values=colors,labels=c("min0maxInf" = 1,"min0max30" = 2, "min2maxInf" = 3,"min2max30" =4, "double" = 5)) +  
  theme(legend.position = 'top') + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Global variant infections\nby day of detection") + theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,0)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_hi) + 
  scale_y_log10(breaks=c(10,100,1000,1e4,1e5,1e6,1e7),labels=trans_format("log10", math_format(10^.x)))+ guides(color=guide_legend(nrow=1))+ theme(legend.key.width = unit(.2,"cm")) + labs(color="Strategy") +
  annotation_logticks(side='l',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) 


### Get global sequencing distributions for each of the strategies
seqrates = cd$Seqrate
seqrates[is.na(seqrates)] = 0

strategy_1_seq_rates = seqrates
strategy_2_seq_rates = strategy_1_seq_rates
strategy_2_seq_rates[strategy_2_seq_rates>30] = 30
strategy_3_seq_rates = seqrates
strategy_3_seq_rates[strategy_3_seq_rates<2] = 2
strategy_4_seq_rates = strategy_3_seq_rates
strategy_4_seq_rates[strategy_4_seq_rates>30] = 30
strategy_5_seq_rates = seqrates*2

popsizes = cd$Population
popsizes[is.na(popsizes)] = 0

strategy_1_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_1_seq_rates[i])))
strategy_2_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_2_seq_rates[i])))
strategy_3_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_3_seq_rates[i])))
strategy_4_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_4_seq_rates[i])))
strategy_5_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_5_seq_rates[i])))

output_df = data.frame(strategy=1:5,tot = c(1,
                                     strategy_2_seq_rates_sum/strategy_1_seq_rates_sum,
                                     strategy_3_seq_rates_sum/strategy_1_seq_rates_sum,
                                     strategy_4_seq_rates_sum/strategy_1_seq_rates_sum,
                                     strategy_5_seq_rates_sum/strategy_1_seq_rates_sum))


# Plot change in time between global detection and local arrival
all_ctrydfs = readRDS("lead_times.rds")
all_ctrydfs = all_ctrydfs %>% group_by(country,strat) %>% dplyr::summarize(mean = mean(mean),continent=continent[1])
all_ctrydfs$strat = factor(all_ctrydfs$strat,levels=c("min0maxInf","min0max30","min2maxInf","min2max30","double"))
all_ctrydfs = all_ctrydfs[all_ctrydfs$country!="Andorra",] # Unrepresentative mobility dynamics
pd = position_dodge(0.5)
lead_by_country_plot = ggplot(all_ctrydfs,aes(x=strat,y=mean)) + 
  theme(legend.key.width = unit(.2,"cm")) + 
  geom_quasirandom(aes(bg=continent,x=strat,y=mean,group=country),pch=21,stroke=0.01,position=pd,cex=.8,alpha=1) + thm +
  ylab("Mean lead time between first global\ndetection and first local case (d)")+ theme(legend.position = 'top')+theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,-20)) +
  scale_fill_npg(labels=c("AF","AS","EU","NA","OC","SA")) + scale_color_npg(labels=c("AF","AS","EU","NA","OC","SA")) + guides(fill=guide_legend(nrow=1),color='none') +
  xlab("") + scale_x_discrete("Strategy",labels=c(1:5)) + labs(fill="Continent") +  theme(legend.title=element_blank()) +
  geom_boxplot(aes(x=strat,y=mean),linewidth=.3,outlier.shape=NA,alpha=0) + guides(bg = guide_legend(override.aes = list(size=1),nrow=1))

# Arrange plot
leg = get_legend(legplot)
ggarrange(
  (time_to_detection_by_strategy|incidence_at_detection_by_strategy|p_local_detection_by_strategy|lead_by_country_plot) + 
    plot_layout(widths=c(0.27,0.27,0.2,0.2)) & plot_annotation(theme = theme(plot.margin = unit(c(0,0,0,0),'cm'))),
  as_ggplot(leg),ncol=1,heights=c(0.9,0.1))

ggsave("figures/Figure_4.pdf",width=180,units="mm",height=45)



####################################
### Specific statistics for text ###
####################################

### Get total expected time to detection and incidence at time of detection
detection_df_grouped = detection_df%>% 
  group_by(strategy) %>% 
  dplyr::summarize(Time_Lo = quantile(detection_day,0.025,na.rm=T),
                   Time_Hi=quantile(detection_day,0.975,na.rm=T),
                   Time_Q25=quantile(detection_day,0.25,na.rm=T),
                   Time_Q75=quantile(detection_day,0.75,na.rm=T),
                   Time_Mean = mean(detection_day,na.rm=T),
                   Inc_Lo = quantile(inc,0.025,na.rm=T),
                   Inc_Hi = quantile(inc,0.975,na.rm=T),
                   Inc_Q25=quantile(inc,0.25,na.rm=T),
                   Inc_Q75=quantile(inc,0.75,na.rm=T),
                   Inc_Mean = mean(inc,na.rm=T),
                   Inc_Median = median(inc,na.rm=T),
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))


### Group outputs by continent for time to detection and compute probability that it is detected in its origin country
detection_df_grouped_continent = detection_df[detection_df$strategy=="min0maxInf",] %>% 
  group_by(onsetcontinent) %>% 
  dplyr::summarize(Time_Lo = quantile(detection_day,0.025,na.rm=T),
                   Time_Hi=quantile(detection_day,0.975,na.rm=T),
                   Time_Q25=quantile(detection_day,0.25,na.rm=T),
                   Time_Q75=quantile(detection_day,0.75,na.rm=T),
                   Time_Mean = mean(detection_day,na.rm=T),
                   Inc_Lo = quantile(inc,0.025,na.rm=T),
                   Inc_Hi = quantile(inc,0.975,na.rm=T),
                   Inc_Mean = mean(inc,na.rm=T),
                   LocalDetect_1=table(continents[detection_country]==onsetcontinent)[1],
                   LocalDetect_0 = table(continents[detection_country]==onsetcontinent)[2],
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))


### Probability of local detection by origin continent
print(data.frame(detection_df_grouped_continent$onsetcontinent,1-detection_df_grouped_continent$LocalDetect_Continent))

### Changes in total global sequencing output
print(1-output_df$tot)


