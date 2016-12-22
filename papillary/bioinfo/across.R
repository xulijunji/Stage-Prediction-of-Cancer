##Contains the code for PCA plots for stages using all genes for all stages as well for individual stages
##using different feature selection technique so far

library(factoextra)
library(FactoMineR)
library(DESeq2)

load('environment/genes_list_DESeq.RData')
load('environment/genes_list_diff_DESeq.RData')
load('environment/dds_tumor_reported.RData')
load('environment/dds.RData')
load('environment/only_tumor_reported.RData')
load('environment/stages_levels.RData')
load('environment/df_stages_rep.RData')
load('environment/genes_list_sam_pair.RData')
load('environment/diff_genes_sam.RData')
load('environment/diff_genes_sam_stages.RData')
load('environment/diff_genes.RData')

combine.all.genes <- function(genes.list)
{
  return(lapply(genes.list, unlist))
}

###DESEQ2
genes.list.comb = combine.all.genes(genes.list)
genes.list.diff.comb = combine.all.genes(genes.list.diff)
genes.list.comb.union = Reduce(union, genes.list.comb)
genes.list.diff.comb.union = Reduce(union, genes.list.diff.comb)
intersect(genes.list.diff.comb.union, genes.list.comb.union)


####SAM
diff.genes.sam.names <- remove.dots(rownames(dds)[diff.genes.sam$`0.05_16`])

genes.stages.sam <- as.numeric(Reduce(union,diff.genes.sam.stages))
genes.stages.sam <- rownames(dds_tumor_reported)[genes.stages.sam]
intersect(genes.stages.sam, genes.list.comb.union)

genes.list.sam.comb <- combine.all.genes(genes.list.sam)
genes.list.sam.comb.union <- Reduce(union, genes.list.sam.comb)
genes.list.sam.comb.union <- rownames(dds_tumor_reported)[genes.list.sam.comb.union]
intersect(genes.list.sam.comb.union, genes.stages.sam)


##Stable Genes
#2nd 150
#1st 200
#cho 200
###################PCA######################################
##Do PCA indiviudal as well as pair wise

pca.comb.union <- PCA(t(only.tumor.reported$dfs$fpqm[genes.list.comb.union,]), graph = F)
fviz_pca_ind(pca.comb.union, habillage = stages.levels, geom = c('point'))
plotPCA(only.tumor.reported$dfs$fpqm[genes.list.diff.comb.union,], 'stage', colData = df.stage.tumor.rep, title = 'comb')

pca.diff.comb.union <- PCA(t(only.tumor.reported$dfs$fpqm[genes.list.diff.comb.union,]), graph = F)
fviz_pca_ind(pca.diff.comb.union, habillage = stages.levels, geom = c('point'))

pca.comb.sam <- PCA(t(only.tumor.reported$dfs$fpqm[genes.stages.sam,]), graph = F)
fviz_pca_ind(pca.comb.sam, habillage = stages.levels, geom = c('point'))

pca.diff.sam <- PCA(t(only.tumor.reported$dfs$fpqm[diff.genes.sam.names,]), graph = F)
fviz_pca_ind(pca.diff.sam, habillage = stages.levels, geom = c('point'))

pca.comb.pair.sam <- PCA(t(only.tumor.reported$dfs$fpqm[genes.list.sam.comb.union,]), graph = F)
fviz_pca_ind(pca.comb.pair.sam, habillage = stages.levels, geom = c('point'))

pca.stable.1 <- PCA(t(only.tumor.reported$dfs$fpqm[stable_genes$fpqm$`1`$V1[1:200],]), graph = F)
fviz_pca_ind(pca.stable.1, habillage = stages.levels, geom = c('point'))

pca.stable.2 <- PCA(t(only.tumor.reported$dfs$fpqm[stable_genes$fpqm$`2`$V1[1:150],]), graph = F)
fviz_pca_ind(pca.stable.2, habillage = stages.levels, geom = c('point'))

pca.stable.cho <- PCA(t(only.tumor.reported$dfs$fpqm[stable_genes$fpqm$`cho`$V1[1:200],]), graph = F)
fviz_pca_ind(pca.stable.cho, habillage = stages.levels, geom = c('point'))
venn.diagram(list(stable_genes$fpqm$`cho`$V1[1:200], stable_genes$fpqm$`1`$V1[1:200],
                  stable_genes$fpqm$`2`$V1[1:150]), category.names = c('cho', '1', '2'),
             filename = 'images/tumor/stable_gene_selection/int.tiff')

####varSelRF genes
pca.all.rf <- PCA(t(only.tumor.reported$dfs$fpqm[all.genes,]), graph = F)
fviz_pca_ind(pca.all.rf, habillage = stages.levels, geom = c('point'))

venn.diagram(list(genes.list.comb.union, genes.list.sam.comb.union, all.genes,
                  stable_genes$fpqm$`2`$V1[1:200]), category.names = c('deseq','sam',
                  'varselrf', 'stable_2'), filename = 'images/tumor/int.tiff')

###stages
genes.varselrf = list()
genes.varselrf[[1]] = list()
genes.varselrf[[1]] = Reduce(union, c(tumor.fpqm.varSelRF.1.2.strat.1$selected.vars,
                                 tumor.fpqm.varSelRF.1.3.strat.1$selected.vars,
                                 tumor.fpqm.varSelRF.1.4$selected.vars))
genes.varselrf[[2]] = Reduce(union, c(tumor.fpqm.varSelRF.1.2.strat.1$selected.vars,
                                      tumor.fpqm.varSelRF.2.3$selected.vars,
                                      tumor.fpqm.varSelRF.2.4$selected.vars))
genes.varselrf[[3]] = Reduce(union, c(tumor.fpqm.varSelRF.2.3.strat.1$selected.vars,
                                      tumor.fpqm.varSelRF.1.3.strat.1$selected.vars,
                                      tumor.fpqm.varSelRF.3.4$selected.vars))
genes.varselrf[[4]] = Reduce(union, c(tumor.fpqm.varSelRF.2.4$selected.vars,
                                      tumor.fpqm.varSelRF.3.4$selected.vars,
                                      tumor.fpqm.varSelRF.1.4$selected.vars))

##Stage i
pca.stage.1.deseq <- PCA(req.dfs$fpqm[stage.ind[[1]],genes.list.comb[[1]]], graph = F)
fviz_pca_ind(pca.stage.1.deseq)
plot(HCPC(pca.stage.1.deseq, 2), choice = 'map')
pca.stage.1.sam <- PCA(req.dfs$fpqm[stage.ind[[1]],genes.list.sam.comb[[1]]], graph = F)
fviz_pca_ind(pca.stage.1.sam)
a = HCPC(pca.stage.1.sam, 2)
plot(a, choice = 'tree')
pca.stage.1.varselrf <- PCA(req.dfs$fpqm[stage.ind[[1]],genes.varselrf[[1]]], graph = F)
fviz_pca_ind(pca.stage.1.varselrf)
HCPC(pca.stage.1.varselrf,2)

#Stage iii
pca.stage.3.deseq <- PCA(req.dfs$fpqm[stage.ind[[3]],genes.list.comb[[3]]], graph = F)
fviz_pca_ind(pca.stage.3.deseq)
HCPC(pca.stage.3.deseq, 3)

pca.stage.3.sam <- PCA(req.dfs$fpqm[stage.ind[[3]],genes.list.sam.comb[[3]]], graph = F)
fviz_pca_ind(pca.stage.3.sam)
HCPC(pca.stage.3.sam, 2)

pca.stage.3.varselrf <- PCA(req.dfs$fpqm[stage.ind[[3]],genes.varselrf[[3]]], graph = F)
fviz_pca_ind(pca.stage.3.varselrf)
HCPC(pca.stage.3.varselrf,2)


####1vs3
pca.1.3.deseq <- PCA(req.dfs$fpqm[Reduce(union,stage.ind[c(1,3)]), genes.list$`stage i`$`stage iii`])
fviz_pca_ind(pca.1.3.deseq, habillage = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
              geom = c('point'))

pca.1.3.sam <- PCA(req.dfs$fpqm[Reduce(union,stage.ind[c(1,3)]), genes.list.sam[[1]]$stageiii])
fviz_pca_ind(pca.1.3.sam, habillage = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
             geom = c('point'), axes = c(1,2))

pca.1.3.varselrf <- PCA(req.dfs$fpqm[Reduce(union,stage.ind[c(1,3)]), tumor.fpqm.varSelRF.1.3.strat.1$selected.vars])
fviz_pca_ind(pca.1.3.varselrf, habillage = droplevels(stages.levels[Reduce(union,stage.ind[c(1,3)])]),
             geom = c('point'), axes = c(1,2))
