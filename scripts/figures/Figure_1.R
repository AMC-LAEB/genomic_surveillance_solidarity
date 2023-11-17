library(ggplot2)
library(patchwork)
library(ggsci)
library(scales)

#Read country information
country_data = readRDS("country_data.rds")
cd = country_data[[2]]
cd = cd[!(is.na(cd$Population)),]

#Plot theme
thm1 =     theme(axis.text = element_text(color="black",size=5),
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
                legend.text=element_text(size=5),
                legend.title =element_text(size=5),
                legend.key=element_rect(fill="white"),
                legend.key.width = unit(.5, "line"),
                legend.spacing.y = unit(.001, 'cm'),
                legend.box.spacing = unit(0, "pt"),
                legend.margin=margin(1,1,1,1))

lty = 1
lw = 0.2
lcol='lightgrey'

# Plot distribution of sequencing rates
sequencing_rate_histogram = ggplot(cd[cd$Seqrate!=0,],aes(x=Seqrate,group=Continent,fill=Continent)) + geom_histogram(linewidth=0.1,col='black',bins=15,alpha=.8) + thm1  +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) + scale_fill_npg(labels=c("AF","AS","EU","NA","OC","SA")) + xlab("Sequencing rate (S/M/wk)") + ylab("Count")+
  theme(legend.position='top') + theme(legend.title=element_blank()) + theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm")) +
  annotation_logticks(side='b',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off')

# Plot distribution of turnaround times
turnaround_time_histogram = ggplot(cd[!(is.na(cd$Median_TAT)),],aes(x=Median_TAT,group=Continent,fill=Continent)) + geom_histogram(linewidth=0.1,col='black',bins=15,alpha=.8) + thm1 +
  scale_fill_npg(labels=c("AF","AS","EU","NA","OC","SA")) + 
  xlab("Median turnaround time (d)") + ylab("Count")  + theme(legend.position='top') + theme(legend.title=element_blank()) +
  theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm"))

# Plot GDP vs sequencing rate
gdp_vs_sequencing_rate = ggplot(cd[!(cd$Seqrate==0) & !(is.na(cd$GDP)),],aes(y=Seqrate,x=GDP)) + geom_point(pch=21,stroke=0.2,aes(group=Continent,fill=Continent),alpha=.8) + thm1 + scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  scale_fill_npg(labels=c("AF","AS","EU","NA","OC","SA")) + scale_x_log10(breaks=c(1000,1e4,1e5),labels=trans_format("log10", math_format(10^.x))) + xlab("GDP per capita (USD)") + ylab("Sequencing rate (S/M/wk)") +
  theme(legend.position='top') + theme(legend.title=element_blank())+ theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm")) + stat_smooth(method='lm',se=F,col='black',lty=2,linewidth=0.4) +
  annotation_logticks(outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off')

# Plot GDP vs median turnaround time
gdp_vs_median_tat = ggplot(cd[!(is.na(cd$Median_TAT)) & !(is.na(cd$GDP)),],aes(y=Median_TAT,x=GDP)) + geom_point(pch=21,stroke=0.2,aes(group=Continent,fill=Continent),alpha=.8) + thm1 + scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  scale_fill_npg(labels=c("AF","AS","EU","NA","OC","SA")) + scale_x_log10(breaks=c(1000,1e4,1e5),labels=trans_format("log10", math_format(10^.x))) + xlab("GDP per capita (USD)") + ylab("Median turnaround time (d)") +
  theme(legend.position='top') + theme(legend.title=element_blank())+ theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm")) + stat_smooth(method='lm',se=F,col='black',lty=2,linewidth=0.2) +
  annotation_logticks(outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off')


# Get cumulative sequencing rates globally
seq_rates_cum = cumsum(sort(unlist(lapply(1:nrow(cd),function(x)rep(cd$Seqrate[x],cd$Population[x]/1e5))),decreasing=TRUE))
inequity_df = data.frame(pop=(1:length(seq_rates_cum))/length(seq_rates_cum),cum = seq_rates_cum/max(seq_rates_cum))

x1_begin = inequity_df$pop[0.5*length(inequity_df$pop)]
x1_end = inequity_df$cum[0.5*length(inequity_df$pop)]
x2_begin = which(inequity_df$cum>0.5)[1]
x2_end = inequity_df$pop[x2_begin]

print("Percentage of population accounting for 50% of sequencing output:")
print(x2_end*100)
print("Percentage of sequencing output accounted for by 50% of the population:")
print(100-x1_end*100)


inequity_df = inequity_df[c(1,seq(1,nrow(inequity_df),100),nrow(inequity_df)),]
inequity_plot = ggplot(inequity_df,aes(x=pop,y=cum)) + geom_line(linewidth=0.4) + thm1 + xlab("Cumulative proportion of global population") + ylab("Cumulative proportion of\nglobal sequencing output") + 
  geom_segment(aes(x=x1_begin,y=0,xend=x1_begin,yend=x1_end),lty=2,linewidth=lw,color=lcol) +
  geom_segment(aes(y=x1_end,x=0,xend=x1_begin,yend=x1_end),lty=2,linewidth=lw,color=lcol) +
  geom_segment(aes(x=0,y=0.5,xend=x2_end,yend=0.5),lty=1,linewidth=lw,color=lcol) +
  geom_segment(aes(x=x2_end,y=0,xend=x2_end,yend=0.5),lty=1,linewidth=lw,color=lcol) +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.008))) + 
  scale_x_continuous(expand=c(0,0.)) + 
  theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm"))
  

# Arrange plots
(sequencing_rate_histogram | inequity_plot) / (turnaround_time_histogram | gdp_vs_sequencing_rate) &
  plot_annotation(theme = theme(plot.margin = unit(c(0,0,0,0),'cm')))
ggsave("figures/Figure_1.pdf",width=90,units="mm",height=80)


gdp_vs_median_tat
ggsave("figures/ExtendedDataFig1.pdf",width=60,units="mm",height=60)


###Flu sequencing rates

# Get cumulative sequencing rates globally
seq_rates_cum = cumsum(sort(unlist(lapply(1:nrow(cd),function(x)rep(cd$Seqrate_Flu[x],cd$Population[x]/1e5))),decreasing=TRUE))
inequity_df = data.frame(pop=(1:length(seq_rates_cum))/length(seq_rates_cum),cum = seq_rates_cum/max(seq_rates_cum))

x1_begin = inequity_df$pop[0.5*length(inequity_df$pop)]
x1_end = inequity_df$cum[0.5*length(inequity_df$pop)]
x2_begin = which(inequity_df$cum>0.5)[1]
x2_end = inequity_df$pop[x2_begin]

inequity_plot_flu = ggplot(inequity_df[c(1,seq(1,nrow(inequity_df),100),nrow(inequity_df)),],aes(x=pop,y=cum)) + geom_line(linewidth=0.4) + thm1 + xlab("Cumulative proportion of global population") + ylab("Cumulative proportion of\nglobal sequencing output") + 
  geom_segment(aes(x=x1_begin,y=0,xend=x1_begin,yend=x1_end),lty=2,linewidth=lw,color=lcol) +
  geom_segment(aes(y=x1_end,x=0,xend=x1_begin,yend=x1_end),lty=2,linewidth=lw,color=lcol) +
  geom_segment(aes(x=0,y=0.5,xend=x2_end,yend=0.5),lty=1,linewidth=lw,color=lcol) +
  geom_segment(aes(x=x2_end,y=0,xend=x2_end,yend=0.5),lty=1,linewidth=lw,color=lcol) +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.008))) + 
  scale_x_continuous(expand=c(0,0.)) + 
  theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm"))


inequity_plot_flu
ggsave("figures/ExtendedDataFig9.pdf",width=60,units="mm",height=60)


## Get correlation coefficients
m = cor.test(log(cd[cd$Seqrate>0,]$Seqrate),log(cd[cd$Seqrate>0,]$GDP),omit.na=T,method='spearman',exact=FALSE)
print(m$estimate)
print(m$p.value)

m = cor.test(log(cd[cd$Seqrate>0,]$Median_TAT),log(cd[cd$Seqrate>0,]$GDP),omit.na=T,method='spearman',exact=FALSE)
print(m$estimate)
print(m$p.value)

# Get number of countries with no sequencing
table(cd$Seqrate==0)

# Get turnaround time quantiles
quantile(cd$Median_TAT,c(0.25,0.75),na.rm=T)
