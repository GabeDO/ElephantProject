library(tidyr)
library(readr)
library(ggplot2)


HetPer10kb <- read_csv("HetPer10kb.csv")

data_long <- gather(HetPer10kb, Elephant, Heterozygosity, 'Female-127_S10':'SRR2009586', factor_key=TRUE)

ggplot(data_long, aes(x=Position, y=Heterozygosity,group=Position))+
  geom_boxplot(notch=FALSE, outlier.shape=NA, fill="red", alpha=0.2)+
  facet_wrap(vars(Chromosome))

ggsave("HetPer10kb_Graph.pdf", width = 25, height = 25)
