library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(corrplot)
library(abind)
library(ggcorrplot)
#load the correlation matrix 

setwd("/data1/2018_ActionDecoding/analysis_fc/OutImage/SeedFC/")
SavePlotDir = "/data1/2018_ActionDecoding/analysis_fc/OutImage/SeedFC/"
#load the matrix, not written here but the location is correct 

#ggplot cannot handle setting the color mapping onto non -1 or 1 so not using ggplot 
# ggcorrplot(OutPlotCor1[indplot,indplot], hc.order = T, type = "upper",
#            outline.col = "white",show.diag = T)
# ggsave(paste(SavePlotDir,"GroupCorrIdentity.eps"),width= 20, height = 20)
# 
# ggcorrplot(OutPlotCor2[indplot,indplot], hc.order = T, type = "upper",
#            outline.col = "white",show.diag = T)
# ggsave(paste(SavePlotDir,"GroupCorrAction.eps"),width= 20, height = 20)
# 
# ggcorrplot(OutPlotCor3[indplot,indplot], hc.order = T, type = "upper",
#            outline.col = "white",show.diag = T)
# ggsave(paste(SavePlotDir,"GroupCorrGoal.eps"),width= 20, height = 20)
# 
# ggcorrplot(OutPlotCor21[indplot,indplot], hc.order = T, type = "upper",
#            outline.col = "white",show.diag = T)
# ggsave(paste(SavePlotDir,"GroupCorrGoal.eps"),width= 20, height = 20)

col4 <- colorRampPalette(c("#00007F","blue" , "#007FFF","cyan",  "#7FFF7F","yellow" , "#FF7F00","red", "#7F0000"))
col1 <- colorRampPalette(c( "#00007F","blue","#007FFF" , "cyan", "white",
                            "yellow", "#FF7F00",  "red","#7F0000"))
col2 <- colorRampPalette(c("blue","white","red"))
pdf(paste(SavePlotDir,"GroupCorrelationPlot_Identity2", ".pdf", sep = ""), width = 20, height = 20)
indplot = c(1:4,7,8)
corrplot(OutPlotCor1[indplot,indplot], diag = T, tl.cex = 1.4, cl.cex = 1.4, method = "color",
         title = "Group Instruction Identity",mar=c(0,0,1,0),
         col=col2(400),cl.lim = c(-1,1),type = "upper",is.corr=T)
dev.off()


pdf(paste(SavePlotDir,"GroupCorrelationPlot_Action2", ".pdf", sep = ""), width = 20, height = 20)
corrplot(OutPlotCor2[indplot,indplot], diag = T, tl.cex = 1.4, cl.cex = 1.4, method = "color",
         title = "Group Instruction Action",mar=c(0,0,1,0),
         col=col2(400),cl.lim = c(-1,1),type = "upper",is.corr=T)
dev.off()
#--ggplot try 


pdf(paste(SavePlotDir,"GroupCorrelationPlot_Goal2", ".pdf", sep = ""), width = 20, height = 20)
corrplot(OutPlotCor3[indplot,indplot], diag = T, tl.cex = 1.4, cl.cex = 1.4, method = "color",
         title = "Group Instruction Goal",mar=c(0,0,1,0),
         col=col2(400),cl.lim = c(-1,1),type = "upper",is.corr=T)
dev.off()



pdf(paste(SavePlotDir,"GroupCorrelationPlot_ActionVSIdentity2", ".pdf", sep = ""), width = 20, height = 20)

# corrplot(OutPlotCor21[1:71,1:71], diag = FALSE, tl.cex = 0.5, method = "color",
#          title = "Group Instruction GoalVSAction",mar=c(0,0,1,0),
#          col=col1(400),is.corr=F)
corrplot(OutPlotCor21[indplot,indplot], diag = T, tl.cex = 1.4, cl.cex = 1.4, method = "color",
         title = "Group Instruction ActionVSIdentity",mar=c(0,0,1,0),
         col=col2(400),cl.lim = c(min(OutPlotCor21[indplot,indplot]),max(OutPlotCor21[indplot,indplot])),type = "lower",is.corr=F)
# ggsave(paste(SavePlotDir,"GroupCorrelationPlot_ActionVSIdentity", ".eps", sep = ""), plot = last_plot(), device = NULL, path = NULL,
#        scale = 1, width = 15, height = 15, units = "in",
#        dpi = 300, limitsize = TRUE)

dev.off()


pdf(paste(SavePlotDir,"GroupCorrelationPlot_GoalVSIdentity2", ".pdf", sep = ""), width = 20, height = 20)
# corrplot(OutPlotCor31[1:71,1:71], diag = FALSE, tl.cex = 0.5, method = "color",
#          title = "Group Instruction GoalVSIdentity",mar=c(0,0,1,0),
#          col=col1(400),is.corr=F)
corrplot(OutPlotCor31[indplot,indplot], diag = T, tl.cex = 1.4,cl.cex = 1.4,  method = "color",
         title = "Group Instruction GoalVSIdentity",mar=c(0,0,1,0),
         col=col2(400),cl.lim =  c(min(OutPlotCor31[indplot,indplot]),max(OutPlotCor31[indplot,indplot])),type = "upper",is.corr=F)
dev.off()


pdf(paste(SavePlotDir,"GroupCorrelationPlot_ActionVSGoal2", ".pdf", sep = ""), width = 20, height = 20)
# corrplot(OutPlotCor21[1:71,1:71], diag = FALSE, tl.cex = 0.5, method = "color",
#          title = "Group Instruction GoalVSAction",mar=c(0,0,1,0),
#          col=col1(400),is.corr=F)
corrplot(OutPlotCor23[indplot,indplot], diag = F, tl.cex = 1.4,cl.cex = 1.4, method = "color",
         title = "Group Instruction GoalVSAction",mar=c(0,0,1,0),
         col=col2(400),cl.lim =  c(min(OutPlotCor23[indplot,indplot]),max(OutPlotCor23[indplot,indplot])),type = "upper",is.corr=F)

dev.off()
