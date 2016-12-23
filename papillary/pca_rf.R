library(randomForest)
library(varSelRF)
library(DESeq2)
library(FactoMineR)
library(MASS)


load('environment/dds_tumor_reported.RData')
load('environment/only_tumor_reported.RData')
load('environment/req_dfs.RData')
load('environment/stages_levels.RData')
load('environment/stages.level.comb.RData')
load('environment/diff_genes.RData')
load('environment/all_genes_varselRF.RData')

library(factoextra)

pca.vs = prcomp(req.dfs$vs)
rf <- randomForest(x = pca.vs$x[,1:5], y = stages.levels)
rf$confusion

pca.vs.diff5 = prcomp(req.dfs$vs[,diff.genes[[5]]])
rf.diff5 <- randomForest(x = pca.vs.diff5$x[,1:10], y = stages.levels)
rf.diff5$confusion

pca.vs.all = prcomp(req.dfs$vs[,all.genes])
rf.all <- randomForest(x = pca.vs.all$x[,1:5], y = stages.levels)
rf.all$confusion

pca.vs.var = prcomp(req.dfs$vs[,only.tumor.reported$genes$vs$var[1:500]])
rf.vs.var <- randomForest(x = pca.vs.var$x[,1:10], y = stages.levels)
rf.vs.var$confusion

rf.comb <- randomForest(x = pca.vs$x[,1:10], y = stages.levels.comb)
rf.comb$confusion

rf.diff5.comb <- randomForest(x = pca.vs.diff5$x[,1:10], y = stages.levels.comb)
rf.diff5.comb$confusion

rf.all.comb <- randomForest(x = pca.vs.all$x[,1:10], y = stages.levels.comb)
rf.all.comb$confusion


lda.dds <- lda(x = t(assay(dds_tumor_reported[diff.genes[[5]],])), grouping = stages.levels
               )
lda.dds.p <- predict(lda.dds)
rf.lda <- randomForest(lda.dds$x, stages.levels)
rf.lda$confusion

train.indexes = c(1:260)
train.indexes = c(1:150,200:240,250:260)
test.indexes = setdiff(c(1:260), train.indexes)

lda.dds.train <- lda(x = req.dfs$vs[train.indexes,all.genes], 
                     grouping = stages.levels.comb[train.indexes])
predictions <- predict(lda.dds.train, newdata = req.dfs$vs[test.indexes,all.genes])
table(stages.levels.comb[test.indexes], predictions$class)

lda.dds.test <- predict(lda.dds.train, 
                        newdata = req.dfs$fpqm[test.indexes,diff.genes[[5]],-train.indexes])
table(lda.dds.test$class)
rf.lda <- randomForest(predictions$x, stages.levels[test.indexes])
rf.lda$confusion

lda.dds <- lda(x = pca.vs$x[train.indexes,1:10], 
                   grouping = droplevels(stages.levels.comb[train.indexes]))
predictions <- predict(lda.dds, newdata = pca.vs$x[test.indexes,1:10])
table(stages.levels.comb[test.indexes], predictions$class)

lda.dds.all <- lda(x = pca.vs.all$x[train.indexes,1:10], 
                   grouping = droplevels(stages.levels[train.indexes]))
predictions <- predict(lda.dds.all, newdata = pca.vs.all$x[test.indexes,1:10])
table(stages.levels[test.indexes], predictions$class)
