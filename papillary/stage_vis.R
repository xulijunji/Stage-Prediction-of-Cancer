#####Inspects the individual stages#########
load('environment/req_dfs.RData')
load('environment/stage.index.RData')

library(DESeq2)
library()
library("RColorBrewer")
library("gplots")
heatmap(req.dfs$vs[stage.ind$`stage ii`,shrunken.genes.vs.1.2.st2])
