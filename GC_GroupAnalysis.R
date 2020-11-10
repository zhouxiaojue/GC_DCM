library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

setwd("/data1/2018_ActionDecoding/analysis_fc/Misc/BME295/")
GC = read.delim('AllSubGCIFGpSTS_RH.txt',header=T)

t.test(GC$IFG2pSTS,GC$pSTS2IFG,paired = T)


ggplot(GC,aes(x=IFG2pSTS,y=pSTS2IFG)) + 
  geom_point()+
  geom_abline(slope =1, intercept = 0) + 
  geom_hline(yintercept = 0.0112,linetype="dotted", color = "blue", size=1.5) +
  geom_vline(xintercept = 0.1363,linetype="dashed", color = "red", size=1.5) +
  labs(title = 'All subjects GC index',x='IFG on pSTS',y='pSTS on IFG')
ggsave("pSTS_RH_IFG_RH_GC.png",width = 4,height = 4)
