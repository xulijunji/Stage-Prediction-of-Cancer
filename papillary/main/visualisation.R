###heatmaps and visualisation based on various features
genes.comb <- list()
gr <- build.groups(length(stages.levels.comb), 5)

genes.group <- pamr.listgene(first.trial$train[[5]], 
                             data = list(x=as.matrix(t(req.dfs$vs[sort(unlist(first.trial$gr[-5])),])),
                                         y=stages.levels.comb[sort(unlist(first.trial$gr[-5]))]), 
                             threshold = first.trial$cv[[5]]$threshold[21],
                             fitcv = first.trial$cv[[5]], genenames = T)
genes.2 <- genes.group[,2][as.numeric(genes.group[,4]) > 0]
genes.1 <- genes.group[,2][as.numeric(genes.group[,3]) > 0]

library(RColorBrewer)
display.brewer.all()
col <- brewer.pal(9, 'PuRd')

create.heatmap(vst_normal_reported[tumor.ind.vs,], stages.levels.comb, 
               get.intersect.genes(g, c(1,2,3,4,5)))
create.heatmap(req.dfs$vs, stages.levels.comb, g2, col = col)
create.heatmap(req.dfs$vs[first.trial$gr[[3]],], 
               stages.levels.comb[first.trial$gr[[3]]], g2)
create.heatmap(vst_normal_reported, 
               colData(dds_tumor_reported_normal_stage)[,4], 
               get.intersect.genes(g, c(1,2,3,4,5)),
               col = col)