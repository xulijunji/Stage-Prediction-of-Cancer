##Generates the genes using indiviudal pair wise results comparisons of stages of DESeq 

install.packages('VennDiagram')
library('VennDiagram')

load('environment/stage_comp.RData')
load('environment/stages_comp_diff.RData')
diff.genes[[1]] = rownames(res)[intersect(which(abs(res$log2FoldChange) > 1),  which(res$padj < 0.01))]
diff.genes[[1]] = remove.dots(diff.genes[[1]])
names(diff.genes) = c(1,2,3,4,5)

stage.comp = list()
stage.comp[['stage i']] <- comp.res(dds_tumor_reported,'stage', 'stage i',c('stage ii', 'stage iii', 'stage iv'))
stage.comp[['stage ii']] <- comp.res(dds_tumor_reported,'stage', 'stage ii',c('stage iii', 'stage iv'))
stage.comp[['stage iii']] <- comp.res(dds_tumor_reported,'stage', 'stage iii',c('stage iv'))

stages.comp.diff <- list()
stages.comp.diff[['stage i']] <- comp.res(dds_tumor_reported[diff.genes[[1]],], 'stage', 'stage i', c('stage ii', 'stage iii', 'stage iv'))
stages.comp.diff[['stage ii']] <- comp.res(dds_tumor_reported[diff.genes[[1]],],'stage', 'stage ii',c('stage iii', 'stage iv'))
stages.comp.diff[['stage iii']] <- comp.res(dds_tumor_reported[diff.genes[[1]],],'stage', 'stage iii',c('stage iv'))



###Using variable from across_tumors.R stage.comp, genes.list and function get.gene
genes.list = list()

##get.genes from across_tumors
################Separating the genes based on original set of genes 
genes.list[['stage i']][['stage ii']] = get.genes(stage.comp$`stage i`$`stage ii`, 2, 0.01, 1)
genes.list[['stage i']][['stage iii']] = get.genes(stage.comp$`stage i`$`stage iii`, 3, 0.01, 1)
genes.list[['stage i']][['stage iv']] = get.genes(stage.comp$`stage i`$`stage iv`, 3, 0.01, 1)
venn.diagram(genes.list$`stage i`, filename = 'bioinfo/stagei/diff2.tiff')
venn.diagram(genes.list$`stage i`, filename = 'bioinfo/stagei/diff2_2.5_2.5.tiff')
venn.diagram(genes.list$`stage i`, filename = 'bioinfo/stagei/diff2_2.5_3.tiff')
venn.diagram(genes.list$`stage i`, filename = 'bioinfo/stagei/diff2_3_3.tiff')
venn.diagram(genes.list$`stage i`, filename = 'bioinfo/stagei/diff2_3_3.5.tiff')
venn.diagram(genes.list$`stage i`, filename = 'bioinfo/stagei/diff2_3.5_3.5.tiff')

genes.list[['stage ii']][['stage i']] = get.genes(stage.comp$`stage i`$`stage ii`, 2, 0.01, 1)
genes.list[['stage ii']][['stage iii']] = get.genes(stage.comp$`stage ii`$`stage iii`, 3, 0.01, 1)
genes.list[['stage ii']][['stage iv']] = get.genes(stage.comp$`stage ii`$`stage iv`, 2.5, 0.01, 1)
venn.diagram(genes.list$`stage ii`, filename = 'bioinfo/stageii/diff2.tiff')
venn.diagram(genes.list$`stage ii`, filename = 'bioinfo/stageii/diff2_2.5_2.5.tiff')
venn.diagram(genes.list$`stage ii`, filename = 'bioinfo/stageii/diff2_2.5_3.tiff')
venn.diagram(genes.list$`stage ii`, filename = 'bioinfo/stageii/diff2_3_2.5.tiff')

genes.list[['stage iii']][['stage i']] = get.genes(stage.comp$`stage i`$`stage iii`, 3.5, 0.01, 1)
genes.list[['stage iii']][['stage ii']] = get.genes(stage.comp$`stage ii`$`stage iii`, 2.5, 0.01, 1)
genes.list[['stage iii']][['stage iv']] = get.genes(stage.comp$`stage iii`$`stage iv`, 2, 0.01, 1)
venn.diagram(genes.list$`stage iii`, filename = 'bioinfo/stageiii/diff2.tiff')
venn.diagram(genes.list$`stage iii`, filename = 'bioinfo/stageiii/diff2.5_2.5_2.tiff')
venn.diagram(genes.list$`stage iii`, filename = 'bioinfo/stageiii/diff3_2.5_2.tiff')
venn.diagram(genes.list$`stage iii`, filename = 'bioinfo/stageiii/diff3_3_2.tiff')
venn.diagram(genes.list$`stage iii`, filename = 'bioinfo/stageiii/diff3.5_2.5_2.tiff')

genes.list[['stage iv']][['stage i']] = get.genes(stage.comp$`stage i`$`stage iv`, 3.0, 0.01, 1)
genes.list[['stage iv']][['stage ii']] = get.genes(stage.comp$`stage ii`$`stage iv`, 3, 0.01, 1)
genes.list[['stage iv']][['stage iii']] = get.genes(stage.comp$`stage iii`$`stage iv`, 2, 0.01, 1)
venn.diagram(genes.list$`stage iv`, filename = 'bioinfo/stageiv/diff2.tiff')
venn.diagram(genes.list$`stage iv`, filename = 'bioinfo/stageiv/diff2.5_2.5_2.tiff')
venn.diagram(genes.list$`stage iv`, filename = 'bioinfo/stageiv/diff3_2.5_2.tiff')
venn.diagram(genes.list$`stage iv`, filename = 'bioinfo/stageiv/diff3_3_2.tiff')


####Separating the genes using only the differentially expressed w.r.t normal and tumor
genes.list.diff <- list()
genes.list.diff[['stage i']][['stage ii']] = get.genes(stages.comp.diff$`stage i`$`stage ii`, 2, 0.01, 1)
genes.list.diff[['stage i']][['stage iii']] = get.genes(stages.comp.diff$`stage i`$`stage iii`, 3, 0.01, 1)
genes.list.diff[['stage i']][['stage iv']] = get.genes(stages.comp.diff$`stage i`$`stage iv`, 3.5, 0.01, 1)
venn.diagram(genes.list.diff$`stage i`, filename = 'images/tumor/DeSeq2/1_fold/stagei/diff2.tiff')
venn.diagram(genes.list.diff$`stage i`, filename = 'images/tumor/DeSeq2/1_fold/stagei/diff2_2.5_2.5.tiff')
venn.diagram(genes.list.diff$`stage i`, filename = 'images/tumor/DeSeq2/1_fold/stagei/diff2_3_3.tiff')
venn.diagram(genes.list.diff$`stage i`, filename = 'images/tumor/DeSeq2/1_fold/stagei/diff2_3_3.5.tiff')

genes.list.diff[['stage ii']][['stage i']] = get.genes(stages.comp.diff$`stage i`$`stage ii`, 2, 0.01, 1)
genes.list.diff[['stage ii']][['stage iii']] = get.genes(stages.comp.diff$`stage ii`$`stage iii`, 3, 0.01, 1)
genes.list.diff[['stage ii']][['stage iv']] = get.genes(stages.comp.diff$`stage ii`$`stage iv`, 3, 0.01, 1)
venn.diagram(genes.list.diff$`stage ii`, filename = 'images/tumor/DeSeq2/1_fold/stageii/diff2.tiff')
venn.diagram(genes.list.diff$`stage ii`, filename = 'images/tumor/DeSeq2/1_fold/stageii/diff2_2.5_2.5.tiff')
venn.diagram(genes.list.diff$`stage ii`, filename = 'images/tumor/DeSeq2/1_fold/stageii/diff2_3_3.tiff')

genes.list.diff[['stage iii']][['stage i']] = get.genes(stages.comp.diff$`stage i`$`stage iii`, 3, 0.01, 1)
genes.list.diff[['stage iii']][['stage ii']] = get.genes(stages.comp.diff$`stage ii`$`stage iii`, 2.5, 0.01, 1)
genes.list.diff[['stage iii']][['stage iv']] = get.genes(stages.comp.diff$`stage iii`$`stage iv`, 2, 0.01, 1)
venn.diagram(genes.list.diff$`stage iii`, filename = 'images/tumor/DeSeq2/1_fold/stageiii/diff2.tiff')
venn.diagram(genes.list.diff$`stage iii`, filename = 'images/tumor/DeSeq2/1_fold/stageiii/diff2_2.5_2.5.tiff')
venn.diagram(genes.list.diff$`stage iii`, filename = 'images/tumor/DeSeq2/1_fold/stageiii/diff3_2.5_3.tiff')

genes.list.diff[['stage iv']][['stage i']] = get.genes(stages.comp.diff$`stage i`$`stage iv`, 3, 0.01, 1)
genes.list.diff[['stage iv']][['stage ii']] = get.genes(stages.comp.diff$`stage ii`$`stage iv`, 2.5, 0.01, 1)
genes.list.diff[['stage iv']][['stage iii']] = get.genes(stages.comp.diff$`stage iii`$`stage iv`, 2, 0.01, 1)
venn.diagram(genes.list.diff$`stage iv`, filename = 'images/tumor/DeSeq2/1_fold/stageiv/diff2.tiff')
venn.diagram(genes.list.diff$`stage iv`, filename = 'images/tumor/DeSeq2/1_fold/stageiv/diff2.5_2.5_2.tiff')
venn.diagram(genes.list.diff$`stage iv`, filename = 'images/tumor/DeSeq2/1_fold/stageiv/diff3_2.5_2.tiff')
