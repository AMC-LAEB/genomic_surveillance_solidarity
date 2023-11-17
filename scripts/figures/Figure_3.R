library(ggridges)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsci)

#Get country data

country_data = readRDS("country_data.rds")
cd = country_data[[1]]
continents = cd$Continent
countries = cd$Country

#Plotting params

width_hi = 1
width_lo = 0.3

fillcol = 'grey95'
nday = 1000


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
            legend.title=element_blank(),
            legend.spacing.x = unit(0.05,"cm"),
            legend.key=element_rect(fill="white"))


### Read genomic surveillance outputs
detection_df = readRDS("detection_outputs.rds")
detection_df = detection_df[detection_df$travel_rate=="mean" & detection_df$phi==1,]

# Keep only scenarios described in the main text
detection_df$scenario = NA
detection_df[detection_df$variant_r0==1.2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.001,]$scenario = 1.2
detection_df[detection_df$variant_r0==1.3 & detection_df$wildtype_beta==0.21 & detection_df$base_prev==0.002,]$scenario = 1.3
detection_df[detection_df$variant_r0==1.6 & detection_df$wildtype_beta==0.22 & detection_df$base_prev==0.005,]$scenario = 1.6
detection_df[detection_df$variant_r0==2 & detection_df$wildtype_beta==0.2 & detection_df$base_prev==0.02,]$scenario = 2

detection_df = detection_df[!(is.na(detection_df$scenario)),]
detection_df = detection_df[!(is.na(detection_df$detection_day)),]

detection_df$onsetcontinent = continents[detection_df$onsetcountry]



### Plot time to detection ridgeplots
time_to_detection = ggplot(detection_df[detection_df$strategy=='min0maxInf' & detection_df$phi == 1,],aes(y=factor(variant_r0),fill=stat(x))) + 
  geom_density_ridges_gradient(rel_min_height=0.0000001,quantile_lines=TRUE,quantiles=c(0.025,0.5,0.975),aes(x=detection_day,group=factor(variant_r0)), 
                               vline_size=0.2, alpha = .3, color = 'black', size=0.1,fill='lightblue') +
  thm + xlab("Day of detection") + ylab(expression("Variant"~ italic("R"[e]))) + xlim(c(0,300)) + theme(plot.margin = unit(c(.2,.1,.2,.1),'cm'))

### Plot incidence by detection ridgeplots
incidence_at_detection = ggplot(detection_df[detection_df$strategy=='min0maxInf' & detection_df$phi == 1,],aes(y=factor(variant_r0),fill=stat(x))) + 
  geom_density_ridges_gradient(rel_min_height=0.001,quantile_lines=TRUE,quantiles=c(0.025,0.5,0.975),aes(x=inc,group=factor(variant_r0)),
                               vline_size=0.2, alpha = .3, color = 'black', size=0.1,fill='lightblue') +thm +
  xlab("Global variant infections\nby day of detection") + ylab(expression("Variant"~ italic("R"[e]))) + scale_x_log10(breaks=c(1,10,100,1000,1e4,1e5,1e6,1e7,1e8),labels=trans_format("log10", math_format(10^.x))) + 
  theme(plot.margin = unit(c(.2,.1,.2,.1),'cm')) + annotation_logticks(side='b',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1)


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
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=.7),linewidth=width_lo) + thm +
  scale_color_npg(labels=c("AF","AS","EU","NA","OC","SA")) + scale_fill_npg(guide='none') + 
  theme(legend.position = 'top') + theme(legend.key.width = unit(.15,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection") + guides(color=guide_legend(nrow=1),fill='none') + 
  theme(legend.margin=margin(-10,0,-10,0),legend.box.margin=margin(-10,-20,-10,-20)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=.7),linewidth=width_hi) +
  theme(plot.margin = unit(c(0,.1,0,.1),'cm'))

### Plot probability of local detection per continent
p_local_detection_per_continent = ggplot(detection_df_grouped[detection_df_grouped$strategy=="min0maxInf",]) + 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) + thm +
  geom_point(aes(x=factor(variant_r0),y=LocalDetect_Continent,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=.7),pch=3,size=1) +
  scale_color_npg(labels=c("AF","AS","EU","NA","OC","SA")) +
  theme(legend.position = 'top')+ theme(legend.key.width = unit(.15,"cm")) + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Day of detection")+ guides(color=guide_legend(nrow=1)) + 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-20,-10,-40)) + ylim(c(0,1))

### Group outputs by continent for incidence at detection
detection_df_grouped = detection_df %>% 
  group_by(onsetcontinent,variant_r0,strategy) %>% 
  dplyr::summarize(Lo = quantile(inc,0.025,na.rm=T),
                   Mean=mean(inc,na.rm=T),
                   Hi=quantile(inc,0.975,na.rm=T),
                   Q25=quantile(inc,0.25,na.rm=T),
                   Q75=quantile(inc,0.75,na.rm=T))

### Plot incidence by detection per continent
incidence_at_detection_per_continent = ggplot(detection_df_grouped[detection_df_grouped$strategy=="min0maxInf",])+ 
  geom_rect(xmin=0,xmax=1.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_rect(xmin=2.5,xmax=3.5,ymin=-20,ymax=300,size=0,fill=fillcol) +
  geom_errorbar(aes(x=factor(variant_r0),ymin=Lo,ymax=Hi,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=.7),linewidth=width_lo) + thm +
  scale_color_npg(labels=c("AF","AS","EU","NA","OC","SA")) + scale_fill_npg(guide='none') + 
  theme(legend.position = 'top') + xlab(expression("Variant"~ italic("R"[e]))) + ylab("Global variant infections\nby day of detection") + theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,-20)) + 
  geom_errorbar(width=0,aes(x=factor(variant_r0),ymin=Q25,ymax=Q75,group=onsetcontinent,col=onsetcontinent),position=position_dodge(width=.7),linewidth=width_hi) + 
  scale_y_log10(breaks=c(1,10,100,1000,1e4,1e5,1e6,1e7),labels=trans_format("log10", math_format(10^.x)))+ guides(color=guide_legend(nrow=1),fill='none')+ theme(legend.key.width = unit(.15,"cm")) + 
  annotation_logticks(side='l',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1)


### Group outputs by onset country
detection_df_grouped_country = detection_df[detection_df$strategy=="min0maxInf",] %>% 
  group_by(onsetcountry,variant_r0,onsetcontinent) %>% 
  dplyr::summarize(Mean = mean(detection_day,na.rm=T),
                   Inc = mean(inc,na.rm=T),
                   LocalDetect_Continent=length(which(continents[detection_country]==onsetcontinent))/length(onsetcontinent))

detection_df_grouped_country$seqrate = cd$Seqrate[detection_df_grouped_country$onsetcountry]
detection_df_grouped_country$seqrate[detection_df_grouped_country$seqrate==0] = NA

### Plot time to detection by origin country's sequencing rate
time_to_detection_by_seqrate = ggplot(detection_df_grouped_country) + 
  labs(color=expression("Variant"~ italic("R"[e]))) + 
  geom_point(aes(x=seqrate,y=Mean,fill=factor(variant_r0)),pch=21,stroke=0.1,alpha=0.3,col='black',cex=.7) + theme(legend.position = 'none')  +
  thm + theme(legend.position='top') + xlab("Origin country sequencing\nrate (S/M/wk)") + ylab("Mean day of detection") +
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) + coord_cartesian(clip = "off") +  scale_size(range=c(0.01,4)) + 
  geom_smooth(aes(x=seqrate,y=Mean,group=factor(variant_r0),color=factor(variant_r0)),method='loess',se=F,linewidth=.5,lty=2) + scale_color_brewer(palette='RdYlBu',direction=-1) + 
  scale_fill_brewer(palette='RdYlBu',direction=-1) + 
  guides(fill='none',col = guide_legend(barheight=.7)) +
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-20,-10,-30),legend.title = element_text(size=5),legend.key.width = unit(.2,"cm")) + ylim(c(0,250)) + 
  annotation_logticks(side='b',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1)

inc_at_detection_by_seqrate = ggplot(detection_df_grouped_country) + 
  labs(color=expression("Variant"~ italic("R"[e]))) + 
  geom_point(aes(x=seqrate,y=Inc,fill=factor(variant_r0)),pch=21,stroke=0.1,alpha=0.3,col='black',cex=.7) + theme(legend.position = 'none')  +
  thm + theme(legend.position='top') + xlab("Origin country sequencing\nrate (S/M/wk)") + ylab("Mean global variant infections\nby day of detection") +
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) + coord_cartesian(clip = "off") +  scale_size(range=c(0.01,4)) + 
  geom_smooth(aes(x=seqrate,y=Inc,group=factor(variant_r0),color=factor(variant_r0)),method='loess',se=F,linewidth=.5,lty=2) + scale_color_brewer(palette='RdYlBu',direction=-1) + 
  scale_fill_brewer(palette='RdYlBu',direction=-1) + 
  guides(fill='none',col = guide_legend(barheight=.7)) + scale_y_log10(breaks=c(1,10,100,1000,1e4,1e5,1e6,1e7),labels=trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(side='bl',outside=F,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + 
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-20,-10,-30),legend.title = element_text(size=5),legend.key.width = unit(.2,"cm"))


### Plot world maps for individual simulations

# Get onset timings for each simulation
getTimings <- function(){

  results = readRDS("simulation_sample.rds")
  nsim = length(results)
  ncountry = dim(results[[1]])[1]
  nday = dim(results[[1]])[2]
  
  onset_days = matrix(NA,nsim,ncountry)
  cuminc = matrix(NA,nsim,nday)
  
  for (i in 1:length(results)){

    sim = results[[i]]
    
    onsets = apply(sim,1,function(x)which(x>=1)[1])
    onset_days[i,] = onsets

  }
  
  return(onset_days)
}

plots_out = list()
sim_idxes = c(74,190,223,208) # Simulations chosen to plot
onset_times = getTimings()

for (q in 1:length(sim_idxes)){
  
  if (q==1){lp='top'} else {lp = 'none'}
  
  plots_out[[q]] = local({
    
    idx= sim_idxes[q]
    d = detection_df[detection_df$variant_r0==1.6 & detection_df$strategy=="min0maxInf",]
    detectiondate = d[idx,]$detection_day
    detectioncountry = d[idx,]$detection_country
    origincountry = d[idx,]$onsetcountry
    countridxes = which(onset_times[idx,]<=detectiondate)
    onsets = onset_times[idx,countridxes]
    
    mp <- NULL
    mapWorld <- borders("world", colour='white', fill="lightgrey", size = .001) # create a layer of borders
    mp <- ggplot() +   mapWorld
    isDetected = countridxes==detectioncountry
    isDetected[countridxes==origincountry] = 2
    lon = cd[countridxes,]$Lon
    lat = cd[countridxes,]$Lat
    
    pl= mp + theme_light()  + 
      theme(axis.line = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank()) + 
      geom_point(aes(size=factor(isDetected),group=factor(isDetected),shape=factor(isDetected),x=as.numeric(lon), y=as.numeric(lat) ,fill=onsets),stroke=0.1) + 
      scale_shape_manual(values=c("0"=21,"1"=24,"2"=25),breaks=c("1","2"),labels=c("Detection country","Origin country")) + 
      scale_size_manual(values=c(1,2,2),guide='none') + 
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust=0.5,size=5),
            axis.ticks = element_blank(),
            legend.title=element_text(size=5),
            legend.key.height = unit(.2,"cm"),
            legend.key.width = unit(.5,"cm"),
            legend.position=lp) + coord_fixed(ylim = c(-50, 85), ratio = 1.3)+ scale_fill_viridis_c(option='C')+ theme(plot.margin=unit(c(-2,-0.1,-1,-0.1),"cm")) + labs(fill="Arrival day") +
      guides(fill=guide_colourbar(title.position="left",title.vjust=0.75,barwidth=3),shape=guide_legend(title.position="left",title.hjust=0.5,title=""))
    print(pl)})
}


### Arrange individual maps
maps =ggarrange(wrap_plots(plots_out) + 
                  plot_layout(guides="collect",ncol=2) & scale_fill_viridis_c(option='C',limits=c(0,100)) 
                & theme(legend.position='bottom',legend.text = element_text(size = 5),legend.box='horizontal',legend.margin=margin(0,0,-10,0),legend.box.margin=margin(-10,-10,-10,-10),legend.spacing.y=unit(.02,'cm')))


### Arrange plots
ggarrange(ggarrange((time_to_detection|incidence_at_detection) + plot_annotation(theme = theme(plot.margin = unit(c(0,.1,0,.1),'cm'))),
                    ((time_to_detection_per_continent|incidence_at_detection_per_continent) + plot_annotation(theme = theme(plot.margin = unit(c(0,.1,0,.1),'cm')))),
                    nrow=1,widths=c(0.38,0.62)),
          ggarrange((time_to_detection_by_seqrate|inc_at_detection_by_seqrate|(p_local_detection_per_continent + ylab("Probability of first detection\nin origin continent"))),
                    maps + plot_annotation(theme = theme(plot.margin = unit(c(-1,.1,0,0),'cm'))),
                    nrow=1,widths=c(0.63,0.37)),ncol=1)

ggsave("figures/Figure_3.pdf",width=180,units="mm",height=80)




