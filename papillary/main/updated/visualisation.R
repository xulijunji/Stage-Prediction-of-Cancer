source('main/updated/initialisation.R')
source('../function.R')
load('environment/sample_info.RData')
load('environment/dds.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')
library(RColorBrewer)

vs_match_tum <- vst(dds)
rownames(vs_match_tum) = remove.dots(rownames(vs_match_tum))
plotPCA(vs_match_tum[diff.genes$`5`,], intgroup = 'type', title = '5 fold genes', colData = colData(dds))
plotPCA(t(vs_normal_comb_reported[, diff.genes$`5`]), intgroup = 'stage.type', title = '5 fold genes', colData = sample.info.all.rep)
plotPCA(t(vst_tumor_tum[,diff.genes$`5`]), intgroup = 'stage.type', title = '5 fold genes', 
        colData = sample.info.all.rep[tumor.ind.vs,])
plotPCA(t(vst_tumor_tum[,net.features.updated$varSelRF$atleast_2]), intgroup = 'stage.type', title = '5 fold genes', 
        colData = data.frame(stage.type=stages.levels.comb))

col <- brewer.pal(9, 'PuRd')
create.heatmap(vs_normal_comb_reported, sample.info.all.rep$stage.type, diff.genes$`5`, 
               'Heatmap using 5 fold genes', col = col, 
               cluster_rows = T, cluster_cols = T)
create.heatmap(vst_tumor_tum, stages.levels.comb, diff.genes$`5`, 
               'Heatmap using 5 fold genes', col = col, 
               cluster_rows = T, cluster_cols = T)
