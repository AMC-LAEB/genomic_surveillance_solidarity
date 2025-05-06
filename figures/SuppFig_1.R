library(ggplot2)
library(patchwork)
library(ggsci)
library(scales)

#Read country information
country_data = readRDS("../outputs/country_data.rds")
cd = country_data[[2]]
cd = cd[!(is.na(cd$Population)),]

#Plot theme
thm1 =     theme(axis.text = element_text(color="black",size=7),
                 panel.spacing.x = unit(1, "mm"),
                 axis.line=element_line(color='black',size = 0.2),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 legend.position='none',
                 axis.title= element_text(size=7),
                 strip.background=element_rect(colour="black",fill="white"),
                 strip.text.x = element_text(size =7),
                 legend.text=element_text(size=7),
                 legend.title =element_text(size=7),
                 legend.key=element_rect(fill="white"),
                 legend.key.width = unit(.5, "line"),
                 legend.spacing.y = unit(.001, 'cm'),
                 legend.box.spacing = unit(0, "pt"),
                 legend.margin=margin(1,1,1,1))

lty = 1
lw = 0.2
lcol='lightgrey'


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

ggsave("SuppFig_1.pdf",width=1000,height=1000,units='px',dpi=320)

