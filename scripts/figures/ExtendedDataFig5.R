library(tidyr)
library(ggplot2)
library(lemon)


# Read analysis output
results = readRDS("GLEAM_validation.rds")

# Plot theme
thm = theme(axis.text = element_text(color="black",size=5),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            #panel.grid.major = element_blank(),
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

### Correlation across all countries
print(cor.test(results[,3],results[,4]))

ggplot(results,aes(x=gleam,y=metapop,group=index_country,fill=index_country)) + geom_point(pch=21,stroke=0.1) + xlab("Median onset (GLEAM)") + ylab("Median onset (metapopulation model)") + thm + 
  geom_abline(lty=2,col='darkgrey') + theme(legend.position='top') + labs(fill="Origin country") + theme(legend.position='none') + facet_rep_wrap(.~index_country,nrow=2,repeat.tick.labels = T)


ggsave("figures/ExtendedDataFig5.pdf",width=180,units="mm",height=90)

