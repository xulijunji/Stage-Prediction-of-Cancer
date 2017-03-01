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
save(train.model, file = 'environment/accuracy_feature/updated/train_model.RData')
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
save(cv.model, file = 'environment/accuracy_feature/updated/cv_model.RData')
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
                                           stages.levels.comb[train.indexes])
save(test.pred, file = 'environment/accuracy_feature/updated/test_pred.RData')   
#############Shrunken features########

###############VarSelRF features###############
train.model[['varSelRF']] <- list()
cv.model[['varSelRF']] <- list()
test.pred[['varSelRF']] <- list()
#############Train########
train.model$varSelRF[['shrunken']] <-build.shrunken.classifier(vst_tumor_tum, train.indexes, 
                                                              net.features.updated$varSelRF[c(3,4,5,6)],
                                                              stages.levels.comb)
train.model$varSelRF[['svm']] <-build.svm.classifier(vst_tumor_tum, train.indexes, 
                                                               net.features.updated$varSelRF[c(3,4,5,6)],
                                                               stages.levels.comb)
train.model$varSelRF[['rf']] <-build.rf.classifier(vst_tumor_tum, train.indexes, 
                                                               net.features.updated$varSelRF[c(3,4,5,6)],
                                                               stages.levels.comb)
train.model$varSelRF[['nb']] <-build.nb.classifier(vst_tumor_tum, train.indexes, 
                                                               net.features.updated$varSelRF[c(3,4,5,6)],
                                                               stages.levels.comb)
########CV########
cv.model$varSelRF[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, net.features.updated$varSelRF[c(3,4,5,6)],
                                               train.model$varSelRF$shrunken, train.indexes, 
                                               stages.levels.comb)
cv.model$varSelRF[['rf']] <- cv.rf.list(vst_tumor_tum, 10, net.features.updated$varSelRF[c(3,4,5,6)],
                                               train.indexes, stages.levels.comb)
cv.model$varSelRF[['svm']] <- cv.svm.list(vst_tumor_tum, 10, net.features.updated$varSelRF[c(3,4,5,6)],
                                        train.indexes, stages.levels.comb)
cv.model$varSelRF[['nb']] <- cv.nb.list(vst_tumor_tum, 10, net.features.updated$varSelRF[c(3,4,5,6)],
                                        train.indexes, stages.levels.comb)
cv.model$varSelRF[['knn']] <- cv.knn.list(vst_tumor_tum, 10, net.features.updated$varSelRF[c(3,4,5,6)],
                                        train.indexes, stages.levels.comb)
########Test#######
test.pred$varSelRF[['shrunken']] <- predict.shrunken(train.model$varSelRF$shrunken, 
                                                     net.features.updated$varSelRF[c(3,4,5,6)],
                                                     vst_tumor_tum[test.indexes, ], 
                                                     cv.model$varSelRF$shrunken$thr)
test.pred$varSelRF[['rf']] <- predict.model(train.model$varSelRF$rf, 
                                            net.features.updated$varSelRF[c(3,4,5,6)],
                                            vst_tumor_tum[test.indexes, ])
test.pred$varSelRF[['svm']] <- predict.model(train.model$varSelRF$svm, 
                                            net.features.updated$varSelRF[c(3,4,5,6)],
                                            vst_tumor_tum[test.indexes, ])
test.pred$varSelRF[['nb']] <- predict.model(train.model$varSelRF$nb, 
                                             net.features.updated$varSelRF[c(3,4,5,6)],
                                             vst_tumor_tum[test.indexes, ])
test.pred$varSelRF[['knn']] <- predict.knn(cv.model$varSelRF$knn, 
                                             net.features.updated$varSelRF[c(3,4,5,6)],
                                             vst_tumor_tum[train.indexes,],
                                             vst_tumor_tum[test.indexes, ], 
                                           stages.levels.comb[train.indexes])
#############VarSelRF features###########

##############DeSeq2 features###########
train.model[[]]