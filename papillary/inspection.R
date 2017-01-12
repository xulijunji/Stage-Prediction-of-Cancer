####A closer look at the four stages####
library(FactoMineR)
library(factoextra)
library(DESeq2)

load('environment/vs_normal_tumor_repored.RData')

genes.most.varying <- rownames(vs.normal.tumor.reported)[get.imp.genes(3,
                                      assay(vs.normal.tumor.reported), 10000)]
sum(rownames(sample.info.all.rep) == colnames(vs.normal.tumor.reported)) == length(rownames(sample.info.all.rep))
indexes.to.remove <- which(sample.info.all.rep$type == 'N')
vs.normal.tumor.reported.varying <- t(assay(vs.normal.tumor.reported)[genes.most.varying,
                                                                      -indexes.to.remove])
sum(rownames(vs.normal.tumor.reported.varying) == sample.info.all.rep$sample.names[sample.info.all.rep$type == 'T'])

rownames(vs.normal.tumor.reported.varying) = c(1:260)
pca.most.var <- PCA(vs.normal.tumor.reported.varying)
fviz_pca_ind(pca.most.var, habillage = stages.levels, 
             geom = c('point'))
res.hcpc <- HCPC(pca.most.var, nb.clust = 4)
