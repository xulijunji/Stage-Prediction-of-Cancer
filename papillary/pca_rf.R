library(randomForest)
library(varSelRF)
library(DESeq2)
library(FactoMineR)
library(MASS)


load('environment/dds_tumor_reported.RData')
load('environment/only_tumor_reported.RData')
load('environment/req_dfs.RData')
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
rf.lda <- randomForest(lda.dds$x, stages.levels)
rf.lda$confusion

train.indexes = c(1:150,200:240,250:260)
test.indexes = setdiff(c(1:260), train.indexes)

lda.dds.train <- lda(x = req.dfs$fpqm[train.indexes,diff.genes[[5]]], 
                     grouping = stages.levels[train.indexes])
predictions <- predict(lda.dds.train, newdata = req.dfs$fpqm[train.indexes,diff.genes[[5]]])
table(stages.levels[train.indexes], predictions$class)

lda.dds.test <- predict(lda.dds.train, 
                        newdata = req.dfs$fpqm[test.indexes,diff.genes[[5]],-train.indexes])
table(lda.dds.test$class)
rf.lda <- randomForest(predictions$x, stages.levels[train.indexes])
rf.lda$confusion
