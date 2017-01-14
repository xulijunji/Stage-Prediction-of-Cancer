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
plot(x = res.hcpc, choice = 'factor', ind.names = T, centers.plot = T,
     new.plot = T)

fit <- hkmeans(vs.normal.tumor.reported.varying, 4)
table(stages.levels, fit$cluster)
fit.kmeans <- kmeans(vs.normal.tumor.reported.varying, 4)
table(stages.levels, fit.kmeans$cluster)

fit.comb.kmeans <- kmeans(req.dfs$vs[,pamr.genes.comb[[5]]], 4)
table(stages.levels.comb, fit.comb.kmeans$cluster)

library(ConsensusClusterPlus)
fit.try <- ConsensusClusterPlus(d = t(vs.normal.tumor.reported.varying),
                                maxK = 4)
View(fit.try[[2]]$consensusMatrix)

