source('main/updated/initialisation.R')

train.model <- list()
cv.model <- list()
test.pred <- list()

##############Shrunken features###############
train.model[['shrunken']] <- list()
cv.model[['shrunken']] <- list()
test.pred[['shrunken']] <- list()
###Shrunken train#########
train.model$shrunken[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, train.indexes, 
                                                          net.features.updated$shrunken[c(3,4,5,6)],
                                                          stages.levels.comb)

train.model$shrunken[['rf']] <- build.rf.classifier(vst_tumor_tum, train.indexes,
                                            net.features.updated$shrunken[c(3,4,5,6)],
                                            stages.levels.comb)

train.model$shrunken[['svm']] <- build.svm.classifier(vst_tumor_tum, train.indexes,
                                                      net.features.updated$shrunken[c(3,4,5,6)],
                                                      stages.levels.comb)

train.model$shrunken[['nb']] <- build.nb.classifier(vst_tumor_tum, train.indexes,
                                                      net.features.updated$shrunken[c(3,4,5,6)],
                                                      stages.levels.comb)
######Shrunken CV######
cv.model$shrunken[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, net.features.updated$shrunken[c(3,4,5,6)],
                                               train.model$shrunken$shrunken, train.indexes, 
                                               stages.levels.comb)
cv.model$shrunken[['svm']] <- cv.svm.list(vst_tumor_tum, 10, net.features.updated$shrunken[c(3,4,5,6)],
                                     train.indexes, stages.levels.comb)
cv.model$shrunken[['rf']] <- cv.rf.list(vst_tumor_tum, 10, net.features.updated$shrunken[c(3,4,5,6)],
                                          train.indexes, stages.levels.comb)
cv.model$shrunken[['nb']] <- cv.nb.list(vst_tumor_tum, 10, net.features.updated$shrunken[c(3,4,5,6)],
                                        train.indexes, stages.levels.comb)
cv.model$shrunken[['knn']] <- cv.knn.list(vst_tumor_tum, 10, net.features.updated$shrunken[c(3,4,5,6)],
                                          train.indexes, stages.levels.comb) 

######Shrunken Test######
test.pred$shrunken[['shrunken']] <- predict.shrunken(train.model$shrunken$shrunken, 
                                                     net.features.updated$shrunken[c(3,4,5,6)],
                                                     vst_tumor_tum[test.indexes, ], 
                                                     cv.model$shrunken$shrunken$thr)
test.pred$shrunken[['rf']] <- predict.model(train.model$shrunken$rf, 
                                                     net.features.updated$shrunken[c(3,4,5,6)],
                                                     vst_tumor_tum[test.indexes, ]) 
test.pred$shrunken[['nb']] <- predict.model(train.model$shrunken$nb, 
                                            net.features.updated$shrunken[c(3,4,5,6)],
                                            vst_tumor_tum[test.indexes, ]) 
test.pred$shrunken[['svm']] <- predict.model(train.model$shrunken$svm, 
                                            net.features.updated$shrunken[c(3,4,5,6)],
                                            vst_tumor_tum[test.indexes, ]) 
test.pred$shrunken[['knn']] <- predict.knn(cv.model$shrunken$knn, 
                                           net.features.updated$shrunken[c(3,4,5,6)],
                                           vst_tumor_tum[train.indexes,], 
                                           vst_tumor_tum[test.indexes, ], 
                                           stages.levels.comb[train.indexes]
                                           )
