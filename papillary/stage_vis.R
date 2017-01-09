#####Inspects the individual stages#########
load('environment/req_dfs.RData')
load('environment/stage.index.RData')

library(DESeq2)
library(pheatmap)
library("RColorBrewer")
library("gplots")
stage.index.1.2 <- Reduce(union, stage.ind[c(1,2)])
pheatmap(req.dfs$vs[stage.index.1.2,shrunken.genes.vs.1.2.st2],
         labels_row = stages.levels[stage.index.1.2],
         )
