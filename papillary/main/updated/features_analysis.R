library(VennDiagram)
load('environment/accuracy_feature/updated/net_features_updated.RData')
load('environment/accuracy_feature/updated/cv_model.RData')
load('environment/accuracy_feature/updated/cv_results.RData')
sapply(net.features.updated$shrunken[c(3,4,5,6)], length)
sapply(net.features.updated$varSelRF[c(3,4,5,6)], length)
lapply(net.features.updated$deseq2[c(3,4,5,6)], function(x)
  {
  sapply(x, length) 
})

create.list.venn <- function(net.features, fold, group)
{
  fold <- as.character(fold)
  group <- as.character(group)
  g1 <- net.features$shrunken[[paste0('atleast_',group)]]
  g2 <- net.features$varSelRF[[paste0('atleast_',group)]]
  g3 <- net.features$deseq2[[paste0('atleast_',group)]][[paste0(fold, ' fold')]]
  req.list <- list(g1,g2,g3)
  names(req.list) <- c('Shrunken', 'VarSelRF', paste0('Deseq2_',fold))
  return(req.list)
}
#1.5 fold
venn.diagram(create.list.venn(net.features.updated, 1.5, 1), filename = 'images/tumor/feature_updated/1_1.5fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1.5, 2), filename = 'images/tumor/feature_updated/2_1.5fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1.5, 3), filename = 'images/tumor/feature_updated/3_1.5fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1.5, 4), filename = 'images/tumor/feature_updated/4_1.5fold.tiff')

#1 fold
venn.diagram(create.list.venn(net.features.updated, 1, 1), filename = 'images/tumor/feature_updated/1_1fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1, 2), filename = 'images/tumor/feature_updated/2_1fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1, 3), filename = 'images/tumor/feature_updated/3_1fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 1, 4), filename = 'images/tumor/feature_updated/4_1fold.tiff')

#2 fold
venn.diagram(create.list.venn(net.features.updated, 2, 1), filename = 'images/tumor/feature_updated/1_2fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 2, 2), filename = 'images/tumor/feature_updated/2_2fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 2, 3), filename = 'images/tumor/feature_updated/3_2fold.tiff')
venn.diagram(create.list.venn(net.features.updated, 2, 4), filename = 'images/tumor/feature_updated/4_2fold.tiff')

####HeatMaps
library(RColorBrewer)
col <- brewer.pal(9, 'PuRd')
create.heatmap(vst_tumor_tum, stages.levels.comb, net.features.updated$deseq2$atleast_4$`2 fold`, 
               '2 fold', col = col)
create.heatmap(vst_tumor_tum, stages.levels.comb, g, 
               '2 fold', col = col)

##Classifiers
temp.class <- list()
temp.class[['shrunken']] <- pamr.train(list(x = t(as.matrix(vst_tumor_tum[train.indexes, g])), 
                                            y = stages.levels.comb[train.indexes]))
temp.class[['svm']] <-  svm(x = vst_tumor_tum[train.indexes, g], y = stages.levels.comb[train.indexes])
temp.class[['nb']] <- naiveBayes(x = vst_tumor_tum[train.indexes, g], y = stages.levels.comb[train.indexes])
temp.class[['rf']] <- randomForest(x = vst_tumor_tum[train.indexes, g], y = stages.levels.comb[train.indexes])



###Predicted
pred.cv.class <- list()
pred.cv.class[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, list(g), list(temp.class$shrunken), 
                                           train.indexes, stages.levels.comb)
pred.cv.class[['knn']] <- find.best.k(vst_tumor_tum[train.indexes, g], 10, stages.levels.comb[train.indexes])

pred.test.class <- list()
pred.test.class[['shrunken']] <- pamr.predict(temp.class$shrunken, t(as.matrix(vst_tumor_tum[test.indexes, g])),
                                     threshold = temp.class$shrunken$threshold[pred.cv.class$shrunken$thr])
pred.test.class[['svm']] <- predict(temp.class$svm, vst_tumor_tum[test.indexes, g])
pred.test.class[['rf']] <- predict(temp.class$rf, vst_tumor_tum[test.indexes, g])
pred.test.class[['nb']] <- predict(temp.class$nb, vst_tumor_tum[test.indexes, g])
pred.test.class[['knn']] <- knn(vst_tumor_tum[train.indexes, g], vst_tumor_tum[test.indexes, g],
                            stages.levels.comb[train.indexes], pred.cv.class$knn$best_k)

##
res.test <- lapply(pred.test.class, function(pred.test){
    get.eval(stages.levels.comb[test.indexes], pred.test)
})

g_1.5_3 <- intersect(intersect(net.features.updated$shrunken$atleast_3, net.features.updated$varSelRF$atleast_3),
               net.features.updated$deseq2$atleast_3$`1.5 fold`) 
g_1.5_2 <- intersect(intersect(net.features.updated$shrunken$atleast_2, net.features.updated$varSelRF$atleast_2),
                     net.features.updated$deseq2$atleast_2$`1.5 fold`) 
g_1.5_1 <- intersect(intersect(net.features.updated$shrunken$atleast_1, net.features.updated$varSelRF$atleast_1),
                     net.features.updated$deseq2$atleast_1$`1.5 fold`) 

g_1_3 <- intersect(intersect(net.features.updated$shrunken$atleast_3, net.features.updated$varSelRF$atleast_3),
                     net.features.updated$deseq2$atleast_3$`1 fold`) 
g_1_2 <- intersect(intersect(net.features.updated$shrunken$atleast_2, net.features.updated$varSelRF$atleast_2),
                   net.features.updated$deseq2$atleast_2$`1 fold`) 
g_1_1 <- intersect(intersect(net.features.updated$shrunken$atleast_1, net.features.updated$varSelRF$atleast_1),
                   net.features.updated$deseq2$atleast_1$`1 fold`) 

g_2_1 <- intersect(intersect(net.features.updated$shrunken$atleast_1, net.features.updated$varSelRF$atleast_1),
                   net.features.updated$deseq2$atleast_1$`2 fold`)
g_2_2 <- intersect(intersect(net.features.updated$shrunken$atleast_2, net.features.updated$varSelRF$atleast_2),
                   net.features.updated$deseq2$atleast_2$`2 fold`)

a <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g, 10)
res_1.5_3 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1.5_3, 10)
res_1.5_2 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1.5_2, 10)
res_1.5_1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1.5_1, 10)

res_1_3 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1_3, 10)
res_1_2 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1_2, 10)
res_1_1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_1_1, 10)

res_2_1 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_2_1, 10)
res_2_2 <- final.res(vst_tumor_tum, train.indexes, test.indexes, stages.levels.comb, g_2_2, 10)

get.conf.mat(a[[2]])
