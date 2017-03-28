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
train.model[['deseq2_2']] <- list()
cv.model[['deseq2_2']] <- list()
test.pred[['deseq2_2']] <- list()

##########Train###############
deseq2_2 <- lapply(net.features.updated$deseq2[c(3,4,5,6)], function(genes.list)
  {
  genes.list[['2 fold']]  
})
deseq2_1.5 <- lapply(net.features.updated$deseq2[c(3,4,5,6)], function(genes.list)
{
  genes.list[['1.5 fold']]  
})
deseq2_1 <- lapply(net.features.updated$deseq2[c(3,4,5,6)], function(genes.list)
{
  genes.list[['1 fold']]  
})
train.model$deseq2_2[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, train.indexes, 
                                                               deseq2_2, stages.levels.comb)
train.model$deseq2_2[['svm']] <- build.svm.classifier(vst_tumor_tum, train.indexes, 
                                                                deseq2_2, stages.levels.comb)
train.model$deseq2_2[['nb']] <- build.nb.classifier(vst_tumor_tum, train.indexes, 
                                                                deseq2_2, stages.levels.comb)
train.model$deseq2_2[['rf']] <- build.rf.classifier(vst_tumor_tum, train.indexes, 
                                                                deseq2_2, stages.levels.comb)
train.model$deseq2_1.5[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, train.indexes, 
                                                                deseq2_1.5, stages.levels.comb)
train.model$deseq2_1.5[['svm']] <- build.svm.classifier(vst_tumor_tum, train.indexes, 
                                                      deseq2_1.5, stages.levels.comb)
train.model$deseq2_1.5[['nb']] <- build.nb.classifier(vst_tumor_tum, train.indexes, 
                                                    deseq2_1.5, stages.levels.comb)
train.model$deseq2_1.5[['rf']] <- build.rf.classifier(vst_tumor_tum, train.indexes, 
                                                    deseq2_1.5, stages.levels.comb)
train.model$deseq2_1[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, train.indexes, 
                                                                deseq2_1, stages.levels.comb)
train.model$deseq2_1[['svm']] <- build.svm.classifier(vst_tumor_tum, train.indexes, 
                                                      deseq2_1, stages.levels.comb)
train.model$deseq2_1[['nb']] <- build.nb.classifier(vst_tumor_tum, train.indexes, 
                                                    deseq2_1, stages.levels.comb)
train.model$deseq2_1[['rf']] <- build.rf.classifier(vst_tumor_tum, train.indexes, 
                                                    deseq2_1, stages.levels.comb)

##########CV###################
cv.model$deseq2_2[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_2,
                                              train.model$deseq2_2$shrunken, train.indexes, 
                                               stages.levels.comb)
cv.model$deseq2_1.5[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_1.5,
                                               train.model$deseq2_1.5$shrunken, train.indexes, 
                                               stages.levels.comb)
cv.model$deseq2_1[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_1,
                                               train.model$deseq2_1$shrunken, train.indexes, 
                                               stages.levels.comb)


cv.model$deseq2_2[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_2,
                                        train.indexes, stages.levels.comb)
cv.model$deseq2_2[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_2,
                                        train.indexes, stages.levels.comb)
cv.model$deseq2_2[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_2,
                                        train.indexes, stages.levels.comb)
cv.model$deseq2_2[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_2,
                                        train.indexes, stages.levels.comb)

cv.model$deseq2_1.5[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_1.5,
                                        train.indexes, stages.levels.comb)
cv.model$deseq2_1.5[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_1.5,
                                        train.indexes, stages.levels.comb)
cv.model$deseq2_1.5[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_1.5,
                                          train.indexes, stages.levels.comb)
cv.model$deseq2_1.5[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_1.5,
                                          train.indexes, stages.levels.comb)

cv.model$deseq2_1[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_1,
                                        train.indexes, stages.levels.comb)
cv.model$deseq2_1[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_1,
                                        train.indexes, stages.levels.comb)
cv.model$deseq2_1[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_1,
                                          train.indexes, stages.levels.comb)
cv.model$deseq2_1[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_1,
                                          train.indexes, stages.levels.comb)


##########Test###############
test.pred$deseq2_2[['shrunken']] <- predict.shrunken(train.model$deseq2_2$shrunken, 
                                                     deseq2_2,
                                                     vst_tumor_tum[test.indexes, ], 
                                                     cv.model$deseq2_2$shrunken$thr)
test.pred$deseq2_2[['svm']] <- predict.model(train.model$deseq2_2$svm, 
                                                     deseq2_2,
                                                     vst_tumor_tum[test.indexes, ])
test.pred$deseq2_2[['rf']] <- predict.model(train.model$deseq2_2$rf, 
                                             deseq2_2,
                                             vst_tumor_tum[test.indexes, ])
test.pred$deseq2_2[['nb']] <- predict.model(train.model$deseq2_2$nb, 
                                             deseq2_2,
                                             vst_tumor_tum[test.indexes, ])
test.pred$deseq2_2[['knn']] <- predict.knn(cv.model$deseq2_2$knn, 
                                    deseq2_2, vst_tumor_tum[train.indexes,],
                                    vst_tumor_tum[test.indexes, ], stages.levels.comb[train.indexes])

test.pred$deseq2_1.5[['shrunken']] <- predict.shrunken(train.model$deseq2_1.5$shrunken, 
                                                     deseq2_1.5,
                                                     vst_tumor_tum[test.indexes, ], 
                                                     cv.model$deseq2_1.5$shrunken$thr)
test.pred$deseq2_1.5[['svm']] <- predict.model(train.model$deseq2_1.5$svm, 
                                             deseq2_1.5,
                                             vst_tumor_tum[test.indexes, ])
test.pred$deseq2_1.5[['rf']] <- predict.model(train.model$deseq2_1.5$rf, 
                                            deseq2_1.5,
                                            vst_tumor_tum[test.indexes, ])
test.pred$deseq2_1.5[['nb']] <- predict.model(train.model$deseq2_1.5$nb, 
                                            deseq2_1.5,
                                            vst_tumor_tum[test.indexes, ])
test.pred$deseq2_1.5[['knn']] <- predict.knn(cv.model$deseq2_1.5$knn, 
                                           deseq2_1.5, vst_tumor_tum[train.indexes,],
                                           vst_tumor_tum[test.indexes, ], stages.levels.comb[train.indexes])

test.pred$deseq2_1[['shrunken']] <- predict.shrunken(train.model$deseq2_1$shrunken, 
                                                     deseq2_1,
                                                     vst_tumor_tum[test.indexes, ], 
                                                     cv.model$deseq2_1$shrunken$thr)
test.pred$deseq2_1[['svm']] <- predict.model(train.model$deseq2_1$svm, 
                                             deseq2_1,
                                             vst_tumor_tum[test.indexes, ])
test.pred$deseq2_1[['rf']] <- predict.model(train.model$deseq2_1$rf, 
                                            deseq2_1,
                                            vst_tumor_tum[test.indexes, ])
test.pred$deseq2_1[['nb']] <- predict.model(train.model$deseq2_1$nb, 
                                            deseq2_1,
                                            vst_tumor_tum[test.indexes, ])
test.pred$deseq2_1[['knn']] <- predict.knn(cv.model$deseq2_1$knn, 
                                           deseq2_1, vst_tumor_tum[train.indexes,],
                                           vst_tumor_tum[test.indexes, ], stages.levels.comb[train.indexes])
save(cv.model, file = 'environment/accuracy_feature/updated/cv_model.RData')
save(test.pred, file = 'environment/accuracy_feature/updated/test_pred.RData')
save(train.model, file = 'environment/accuracy_feature/updated/train_model.RData')
