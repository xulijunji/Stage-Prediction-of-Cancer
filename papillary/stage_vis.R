#####Inspects the individual stages#########
load('environment/req_dfs.RData')
load('environment/stage.index.RData')

library(DESeq2)
heatmap(req.dfs$vs[stage.ind$`stage ii`, 1:10])
