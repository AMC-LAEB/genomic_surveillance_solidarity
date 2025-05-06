library(ggplot2)
library(ggsci)
library(latex2exp)
library(lemon)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)

# Get country data

country_data = readRDS("../outputs/country_data.rds")
cd = country_data[[1]]
continents = cd$Continent
countries = cd$Country


# Plotting parameters

thm = theme(axis.text = element_text(color="black",size=7),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position='none',
            axis.title= element_text(size=7),
            strip.background=element_rect(colour="black",fill="white"),
            strip.text.x = element_text(size = 5),
            legend.text = element_text(size = 7),
            legend.title=element_blank(),
            legend.spacing.x = unit(0.05,"cm"),
            legend.key=element_rect(fill="white"))

fillcol = 'grey90'
nday = 1000

width_hi = 1
width_lo = 0.3

f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
cols = f('Set1')[c(1:5,7)]
cols_strategies = rev(c("#7B3F00","navy","darkgreen","red",'#008080'))

### Read genomic surveillance outputs
#detection_df = readRDS("detection_outputs_threshold_F.rds")
detection_df = readRDS("../outputs/combined_outputs.rds")
detection_df = detection_df[detection_df$strategy!="double",]
detection_df[detection_df$detection_day==0,]$detection_day = NA


detection_df = detection_df[detection_df$travel_rate=="mean",]
detection_df = detection_df[detection_df$phi == 1,]

detection_df$strategy = factor(detection_df$strategy,levels=c("min0maxInf","min2maxInf","min0max30","min2max30"))


# Keep only scenarios described in the main text
detection_df$scenario = NA
detection_df[detection_df$variant_r0==1.2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.001,]$scenario = 1.2
detection_df[detection_df$variant_r0==1.3 & detection_df$wildtype_beta==0.21 & detection_df$base_prev==0.002,]$scenario = 1.3
detection_df[detection_df$variant_r0==1.6 & detection_df$wildtype_beta==0.22 & detection_df$base_prev==0.005,]$scenario = 1.6
detection_df[detection_df$variant_r0==2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.02,]$scenario = 2
detection_df = detection_df[!(is.na(detection_df$scenario)),]
detection_df = detection_df[!(is.na(detection_df$detection_day)),]

detection_df$onsetcontinent = continents[detection_df$onset_country]
detection_df$detectioncontinent = continents[detection_df$detection_country]


### Group time to detection by r0 and surveillance strategy
detection_df_grouped = detection_df %>% group_by(variant_r0,strategy) %>% dplyr::summarize(Lo = quantile(detection_day,0.025,na.rm=T),
                                                                                           Hi=quantile(detection_day,0.975,na.rm=T),
                                                                                           Mean = mean(detection_day,na.rm=T),
                                                                                           Q25=quantile(detection_day,0.25,na.rm=T),
                                                                                           Q75=quantile(detection_day,0.75,na.rm=T),
                                                                                           LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))

detection_df_grouped$strategy = factor(detection_df_grouped$strategy,levels=c("min0maxInf","min2maxInf","min0max30","min2max30"))

### Plot incidence at detection by strategy by r0
time_to_detection_by_strategy = ggplot(detection_df_grouped) + 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=330,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=330,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=0.85),linewidth=width_lo) + thm +
  #scale_color_manual(values=cols_strategies) + scale_fill_aaas(palette = "Dark2",guide='none') + 
  theme(legend.position = 'top') + theme(legend.key.width = unit(.2,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection") + 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(10,0,10,0)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=strategy,col=strategy),position=position_dodge(width=0.85),linewidth=width_hi) +
  scale_color_manual(values=cols_strategies,labels=c("2022","A","B","C","D")) + 
  scale_fill_manual(values=cols_strategies) +
  geom_point(aes(x=factor(variant_r0),y=Mean,group=strategy,col=strategy),position=position_dodge(width=0.85),cex=1,pch=21,fill='white') +
  guides(color=guide_legend(nrow=1,override.aes = list(shape = NA)),fill='none') + theme(legend.key.size=unit(width_hi,'lines')) + labs(color="Strategy") +theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,0))


### Get legend
legplot = ggplot(detection_df_grouped)+ 
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=.7),linewidth=width_hi) + thm +
  scale_color_manual(values=cols_strategies) + scale_fill_manual(values=cols_strategies,guide='none') + 
  scale_color_manual(values=cols_strategies,labels=c("2022: 2022 baseline","A: 2022 baseline + minimum global capacity",  "B: 2022 baseline, capped at 30 S/M/wk","C: 2022 baseline, capped at 30 S/M/wk\n+ minimum global capacity")) +
  scale_fill_aaas() +
  guides(color=guide_legend(nrow=3),fill='none') + theme(legend.key.size=unit(.6,'lines')) + labs(color="Strategy") + theme(legend.position = 'top')+ theme(legend.key.width = unit(.2,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection")+ 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,0))

### Plot probability of local detection by strategy by r0
p_local_detection_by_strategy = ggplot(detection_df_grouped)+ 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_point(aes(x=factor(variant_r0),y=LocalDetect_Continent,group=strategy,col=strategy),position=position_dodge(width=.7),pch=3,size=1) + thm +
  scale_color_manual(values=cols_strategies,labels=c("min0maxInf" = "2022","min2maxInf" = "A","min0max30" = "B", "min2max30" = "C")) + scale_fill_aaas(guide='none') + 
  theme(legend.position = 'top')+ theme(legend.key.width = unit(.2,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Probability of first detection\nin origin continent")+ 
  guides(color=guide_legend(nrow=1)) + 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,10,-10,0)) + labs(color="Strategy") + ylim(c(0.6,1)) 


### Group incidence by r0 and strategy
detection_df_grouped = detection_df %>% group_by(variant_r0,strategy) %>% dplyr::summarize(Lo = quantile(detection_infections,0.025,na.rm=T),
                                                                                           Hi=quantile(detection_infections,0.975,na.rm=T),
                                                                                           Med = quantile(detection_infections,0.5,na.rm=T),
                                                                                           Q25=quantile(detection_infections,0.25,na.rm=T),
                                                                                           Q75=quantile(detection_infections,0.75,na.rm=T),
                                                                                           Mean=mean(detection_infections,na.rm=T))

# Plot incidence at detection by strategy
detection_df_grouped$strategy = factor(detection_df_grouped$strategy,levels=c("min0maxInf","min2maxInf","min0max30","min2max30"))
incidence_at_detection_by_strategy = ggplot(detection_df_grouped) + 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=strategy,col=strategy),position=position_dodge(width=0.85),linewidth=width_lo) + thm +
  scale_color_manual(values=cols_strategies,labels=c("min0maxInf" = "2022","min2maxInf" = "A", "min0max30" = "B","min2max30" ="C")) + scale_fill_brewer(palette = "Accent",guide='none') + 
  scale_fill_manual(values=cols_strategies,guide='none') + 
  theme(legend.position = 'top') + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Variant infections by detection day") + theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,0)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=strategy,col=strategy),position=position_dodge(width=0.85),linewidth=width_hi) + 
  scale_y_log10(breaks=c(10,100,1000,1e4,1e5,1e6,1e7),labels=trans_format("log10", math_format(10^.x)))+ guides(color=guide_legend(nrow=1,override.aes = list(shape = NA)))+ theme(legend.key.width = unit(.2,"cm")) + labs(color="Strategy") +
  annotation_logticks(side='l',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) +
  geom_point(aes(x=factor(variant_r0),y=Mean,group=strategy,col=strategy),position=position_dodge(width=0.85),cex=1,pch=21,fill='white') 



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

popsizes = cd$Population
popsizes[is.na(popsizes)] = 0

strategy_1_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_1_seq_rates[i])))
strategy_2_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_2_seq_rates[i])))
strategy_3_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_3_seq_rates[i])))
strategy_4_seq_rates_sum = sum(unlist(lapply(1:length(cd$Seqrate),function(i)cd$Population[i]/1e6* strategy_4_seq_rates[i])))

output_df = data.frame(strategy=1:4,tot = c(1,
                                            strategy_2_seq_rates_sum/strategy_1_seq_rates_sum,
                                            strategy_3_seq_rates_sum/strategy_1_seq_rates_sum,
                                            strategy_4_seq_rates_sum/strategy_1_seq_rates_sum))

# Plot global sequencing output by strategy                             
distribution_by_strategy =  ggplot(output_df,aes(x=factor(strategy),fill=factor(strategy),y=tot)) + geom_bar(stat='identity',width=0.8) + thm + 
  scale_fill_aaas(labels=c("2022","A","B","C")) + xlab("Strategy") + ylab("Global sequencing output\n relative to 2022 baseline") + 
  scale_y_continuous(breaks=seq(0,1,0.25))


# Plot change in time between global detection and local arrival
all_ctrydfs = readRDS("../outputs/lead_times.rds")
all_ctrydfs = all_ctrydfs[all_ctrydfs$strat!="double",]

all_ctrydfs[all_ctrydfs$country=="United States of America" & all_ctrydfs$r0==1.6,]
all_ctrydfs[all_ctrydfs$country%in% c("Rwanda","Kazakhstan","Indonesia","United Kingdom") & 
              all_ctrydfs$r0==1.3 & all_ctrydfs$strat %in% c("min0maxInf","min2max30"),]

all_ctrydfs = all_ctrydfs %>% group_by(country,strat) %>% dplyr::summarize(mean = mean(mean),continent=continent[1])
all_ctrydfs$strat = factor(all_ctrydfs$strat,levels=c("min0maxInf","min2maxInf","min0max30","min2max30"))
all_ctrydfs = all_ctrydfs[all_ctrydfs$country!="Andorra",] # Unrepresentative mobility dynamics
pd = position_dodge(0.5)

lead_by_country_plot = ggplot(all_ctrydfs,aes(x=strat,y=mean)) + 
  theme(legend.key.width = unit(.2,"cm")) + 
  geom_quasirandom(aes(bg=continent,x=strat,y=mean,group=country),pch=21,stroke=0.01,position=pd,cex=.8,alpha=1) + thm +
  ylab("Mean lead time (d)")+ theme(legend.position = 'top')+theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,-20)) +
  scale_fill_manual(values=cols,labels=c("AF","AS","EU","NA","OC","SA")) + scale_color_manual(values=cols,labels=c("AF","AS","EU","NA","OC","SA")) + guides(fill=guide_legend(nrow=1),color='none') +
  xlab("") + scale_x_discrete("Strategy",labels=c("2022","A","B","C")) + labs(fill="Continent") +  theme(legend.title=element_blank()) +
  geom_boxplot(aes(x=strat,y=mean),linewidth=.3,outlier.shape=NA,alpha=0) + guides(bg = guide_legend(override.aes = list(size=1.5),nrow=1))

# Arrange plot
leg = get_legend(legplot+theme(legend.title=element_text(size=7)))

time_to_detection_by_strategy+theme(plot.margin=unit(c(0,0.1,0,-0.1),"cm"))+
  incidence_at_detection_by_strategy+theme(plot.margin=unit(c(0,0.1,0,-0.1),"cm"))+
  p_local_detection_by_strategy+theme(plot.margin=unit(c(0,0.1,0,-0.1),"cm"))+
  lead_by_country_plot + theme(plot.margin=unit(c(0,0.15,0,-0.1),"cm"))+
  plot_layout(ncol=4) #& plot_annotation(theme = theme(plot.margin = unit(c(0,0,0,0),'cm')))
ggsave("Figure_3.pdf",width=2080,height=650,units='px',dpi=320)



#################################### 
### Specific statistics for text ###
####################################

### Get total expected time to detection and incidence at time of detection
detection_df_grouped = detection_df%>% 
  group_by(strategy) %>% 
  dplyr::summarize(Time_Mean = mean(detection_day,na.rm=T),
                   Time_Lo = quantile(detection_day,0.025,na.rm=T),
                   Time_Hi=quantile(detection_day,0.975,na.rm=T),
                   Time_Q25=quantile(detection_day,0.25,na.rm=T),
                   Time_Q75=quantile(detection_day,0.75,na.rm=T),
                   Inc_Mean = mean(detection_infections,na.rm=T),
                   Inc_Lo = quantile(detection_infections,0.025,na.rm=T),
                   Inc_Hi = quantile(detection_infections,0.975,na.rm=T),
                   Inc_Q25=quantile(detection_infections,0.25,na.rm=T),
                   Inc_Q75=quantile(detection_infections,0.75,na.rm=T),
                   Inc_Median = median(detection_infections,na.rm=T),
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))


### Group outputs by continent for time to detection and compute probability that it is detected in its origin country
detection_df_grouped_continent = detection_df[detection_df$strategy=="min0maxInf",] %>% 
  group_by(onsetcontinent) %>% 
  dplyr::summarize(Time_Lo = quantile(detection_day,0.025,na.rm=T),
                   Time_Hi=quantile(detection_day,0.975,na.rm=T),
                   Time_Q25=quantile(detection_day,0.25,na.rm=T),
                   Time_Q75=quantile(detection_day,0.75,na.rm=T),
                   Time_Mean = mean(detection_day,na.rm=T),
                   Inc_Lo = quantile(detection_infections,0.025,na.rm=T),
                   Inc_Hi = quantile(detection_infections,0.975,na.rm=T),
                   Inc_Mean = mean(detection_infections,na.rm=T),
                   LocalDetect_1=table(continents[detection_country]==onsetcontinent)[1],
                   LocalDetect_0 = table(continents[detection_country]==onsetcontinent)[2],
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))


### Probability of local detection by origin continent
print(data.frame(detection_df_grouped_continent$onsetcontinent,1-detection_df_grouped_continent$LocalDetect_Continent))

### Changes in total global sequencing output
print(1-output_df$tot)


