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
