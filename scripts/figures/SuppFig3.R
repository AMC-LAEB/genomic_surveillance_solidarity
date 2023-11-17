country_data = readRDS("country_data.rds") 
travel_rates = colSums(country_data[[4]],na.rm=T)
plot(log(travel_rates),log(country_data[[1]]$Seqrate))
seqrate = country_data[[1]]$Seqrate

cor.test(log(travel_rates[seqrate>0]),log(seqrate[seqrate>0]),method='spearman',exact=F)$p.value
cor.test(log(travel_rates[seqrate>0]),log(seqrate[seqrate>0]),method='spearman',exact=F)$estimate



df = data.frame(seqrate = seqrate[seqrate>0], travelrate = travel_rates[seqrate>0],continent=countrydata[[1]]$Continent[seqrate>0])

df = df[df$travelrate>1e-7,]

ggplot(df,aes(y=seqrate,x=travelrate)) + geom_point(pch=21,stroke=0.2,aes(group=continent,fill=continent),alpha=0.8) + thm1 + scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),labels=trans_format("log10", math_format(10^.x))) +
  scale_fill_npg(labels=c("AF","AS","EU","NA","OC","SA")) + scale_x_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1),labels=trans_format("log10", math_format(10^.x))) + xlab("Daily per capita mobility rate") + ylab("Sequencing rate (S/M/wk)") +
  theme(legend.position='top') + theme(legend.title=element_blank())+ theme(legend.key.width=unit(.2,"cm")) + theme(legend.key.height=unit(.2,"cm")) + stat_smooth(method='lm',se=F,col='black',lty=2,linewidth=0.4) +
  annotation_logticks(outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off')

ggsave("figures/SuppFig3.pdf",width=60,units="mm",height=60)
