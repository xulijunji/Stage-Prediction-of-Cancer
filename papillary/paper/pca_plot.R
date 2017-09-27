setwd('~/h_short/papillary/')
source('../function.R')
source('../Feature_Extract.R')
load('environment/dds.RData')
load('environment/accuracy_feature/vs_normal_tumor_repored.RData')
load('environment/sample_info.RData')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/diff_genes.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')

library(DESeq2)
rownames(dds) <- remove.dots(rownames(dds))
p1 <- plotPCA(vst(dds)[diff.genes$`5`,], intgroup = 'type', title = 'Matched normal and tumor samples')
p2 <- plotPCA(t(vs_normal_comb_reported[,diff.genes$`5`]), intgroup = 'stage.type', title = 'Normal and all tumor samples with reported stages',
              colData = sample.info.all.rep)
multiplot(list(p1, p2), cols = 2)

load('environment/genes_list_DESeq.RData')
load('environment/genes_list_diff_DESeq.RData')
load('environment/accuracy_feature/tumor_ind_vs.RData')
load('environment/stages.level.comb.RData')
g <- unique(unlist(g))
p3 <- plotPCA(t(vs_normal_comb_reported[tumor.ind.vs,g]), intgroup = 'stage.type',
              title = '3A', colData = sample.info.all.rep[tumor.ind.vs,]
                             )
load('environment/stage_comp.RData')
g <- list()
g[['stage i']] <- list()
g[['stage i']][['stage ii']] = get.genes(stage.comp$`stage i`$`stage ii`, 2, 0.01, 1)
g[['stage i']][['stage iii']] = get.genes(stage.comp$`stage i`$`stage iii`, 2, 0.01, 1)
g[['stage i']][['stage iv']] = get.genes(stage.comp$`stage i`$`stage iv`, 2, 0.01, 1)

g[['stage ii']] <- list()
g[['stage ii']][['stage iii']] = get.genes(stage.comp$`stage ii`$`stage iii`, 2, 0.01, 1)
g[['stage ii']][['stage iv']] = get.genes(stage.comp$`stage ii`$`stage iv`, 2, 0.01, 1)

g[['stage iii']] <- list()
g[['stage iii']][['stage iv']] = get.genes(stage.comp$`stage iii`$`stage iv`, 2, 0.01, 1)

levels(stages.levels.comb) <- c('Early Stage', 'Late Stage')
d <- data.frame(stage=stages.levels.comb, row.names = colnames(vs.normal.tumor.reported[,tumor.ind.vs]))
p4 <- plotPCA(t(vs_normal_comb_reported[tumor.ind.vs, g]), intgroup = 'stage',
              title = '3B', colData = d)
multiplot(list(p3, p4), col = 2)              
