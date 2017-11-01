setwd('~/h_short/papillary/')
source('../function.R')
source('../Feature_Extract.R')
source('paper/multiplot.R')
load('environment/dds.RData')
load('environment/accuracy_feature/vs_normal_tumor_repored.RData')
load('environment/sample_info.RData')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/diff_genes.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')
load('environment/accuracy_feature/updated/vst_tum_rep.RData')

library(DESeq2)
rownames(dds) <- remove.dots(rownames(dds))
p1 <- plotPCA(vst(dds)[diff.genes$`5`,], intgroup = 'type', title = 'Matched normal and tumor samples')
s1.list <- lapply(c(1,2,3,4), function(i)
  {
  plotPCA(vst(dds)[diff.genes[[i]],], intgroup = 'type',
          title=bquote(log[2]~FC~.(names(diff.genes[i]))))
          #title = bquote(a~[2]~FC~.(names(diff.genes)[i]))
})
multiplot(s1.list, cols = 2)
s2.a <- plotPCA(vst(dds)[diff.genes.up$`5`,], intgroup = 'type',
                title = 'A')
s2.b <- plotPCA(vst(dds)[diff.genes.low$`5`,], intgroup = 'type',
                title = 'B')
s2.b

grid_arrange_shared_legend(list(s2.a, s2.b))
p2 <- plotPCA(t(vs_normal_comb_reported[,diff.genes$`5`]), intgroup = 'stage.type', title = 'Normal and all tumor samples with reported stages',
              colData = sample.info.all.rep)

p1.up <- plotPCA(vst(dds)[diff.genes.up$`2`,], intgroup = 'type', title = 'Matched normal and tumor samples')
multiplot(list(p1, p2), cols = 2)

load('environment/genes_list_DESeq.RData')
load('environment/genes_list_diff_DESeq.RData')
load('environment/accuracy_feature/tumor_ind_vs.RData')
load('environment/stages.level.comb.RData')
g <- unique(unlist(g))
p3 <- plotPCA(t(vs_normal_comb_reported[tumor.ind.vs,g]), intgroup = 'stage.type',
              title = '3A', colData = sample.info.all.rep[tumor.ind.vs,])

                                           
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
