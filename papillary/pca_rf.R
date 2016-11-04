library(randomForest)
library(varSelRF)
library(DESeq2)
library(FactoMineR)
library(MASS)


load('environment/dds_tumor_reported.RData')
load('environment/stages_levels.RData')
load('environment/diff_genes.RData')

pca.dds = prcomp(t(assay(dds_tumor_reported)))
rf <- randomForest(x = pca.dds$x[,1:10], y = stages.levels)
rf$confusion

pca.dds.diff5 = prcomp(t(assay(dds_tumor_reported[diff.genes[[5]],])))
rf.diff5 <- randomForest(x = pca.dds.rem$x[,1:10], y = stages.levels)
rf.diff$confusion


lda.dds <- lda(x = t(assay(dds_tumor_reported[diff.genes[[5]],])), grouping = stages.levels
               )
lda.dds.p <- predict(lda.dds)
rf.lda <- randomForest(lda.dds.p$x, stages.levels)
rf.lda$confusion

train.indexes = c(c(1:150,200:240,250:260))
lda.dds.train <- lda(x = t(assay(dds_tumor_reported[diff.genes[[5]],train.indexes])), 
                     grouping = stages.levels[train.indexes])
lda.dds.test <- predict(lda.dds.train, 
                        newdata = t(assay(dds_tumor_reported[diff.genes[[5]],])))
rf.lda <- randomForest(lda.dds.test$x, stages.levels)
rf.lda$confusion
