source('main/updated/initialisation.R')

train.model <- list()
cv.model <- list()
test.pred <- list()

##############Shrunken features###############
tr.ind <- train.trial.ind
te.ind <- test.trial.ind
gr.final <- gr.trial
net.fea <- net.features.trial

train.model[['shrunken']] <- list()
cv.model[['shrunken']] <- list()
test.pred[['shrunken']] <- list()
###Shrunken train#########
train.model$shrunken[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                          net.fea$shrunken[c(3,4,5,6)],
                                                          stages.levels.comb)

train.model$shrunken[['rf']] <- build.rf.classifier(vst_tumor_tum, tr.ind,
                                            net.fea$shrunken[c(3,4,5,6)],
                                            stages.levels.comb)

train.model$shrunken[['svm']] <- build.svm.classifier(vst_tumor_tum, tr.ind,
                                                      net.fea$shrunken[c(3,4,5,6)],
                                                      stages.levels.comb)

train.model$shrunken[['nb']] <- build.nb.classifier(vst_tumor_tum, tr.ind,
                                                      net.fea$shrunken[c(3,4,5,6)],
                                                      stages.levels.comb)
save(train.model, file = 'environment/accuracy_feature/updated/train_model.RData')
######Shrunken CV######
cv.model$shrunken[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, net.fea$shrunken[c(3,4,5,6)],
                                               train.model$shrunken$shrunken, tr.ind, 
                                               stages.levels.comb)
cv.model$shrunken[['svm']] <- cv.svm.list(vst_tumor_tum, 10, net.fea$shrunken[c(3,4,5,6)],
                                     tr.ind, stages.levels.comb)
cv.model$shrunken[['rf']] <- cv.rf.list(vst_tumor_tum, 10, net.fea$shrunken[c(3,4,5,6)],
                                          tr.ind, stages.levels.comb)
cv.model$shrunken[['nb']] <- cv.nb.list(vst_tumor_tum, 10, net.fea$shrunken[c(3,4,5,6)],
                                        tr.ind, stages.levels.comb)
cv.model$shrunken[['knn']] <- cv.knn.list(vst_tumor_tum, 10, net.fea$shrunken[c(3,4,5,6)],
                                          tr.ind, stages.levels.comb) 
save(cv.model, file = 'environment/accuracy_feature/updated/cv_model.RData')
######Shrunken Test######
test.pred$shrunken[['shrunken']] <- predict.shrunken(train.model$shrunken$shrunken, 
                                                     net.fea$shrunken[c(3,4,5,6)],
                                                     vst_tumor_tum[te.ind, ], 
                                                     cv.model$shrunken$shrunken$thr)
test.pred$shrunken[['rf']] <- predict.model(train.model$shrunken$rf, 
                                                     net.fea$shrunken[c(3,4,5,6)],
                                                     vst_tumor_tum[te.ind, ]) 
test.pred$shrunken[['nb']] <- predict.model(train.model$shrunken$nb, 
                                            net.fea$shrunken[c(3,4,5,6)],
                                            vst_tumor_tum[te.ind, ]) 
test.pred$shrunken[['svm']] <- predict.model(train.model$shrunken$svm, 
                                            net.fea$shrunken[c(3,4,5,6)],
                                            vst_tumor_tum[te.ind, ]) 
test.pred$shrunken[['knn']] <- predict.knn(cv.model$shrunken$knn, 
                                           net.fea$shrunken[c(3,4,5,6)],
                                           vst_tumor_tum[tr.ind,], 
                                           vst_tumor_tum[te.ind, ], 
                                           stages.levels.comb[tr.ind])
save(test.pred, file = 'environment/accuracy_feature/updated/test_pred.RData')   
#############Shrunken features########

###############VarSelRF features###############
train.model[['varSelRF']] <- list()
cv.model[['varSelRF']] <- list()
test.pred[['varSelRF']] <- list()
#############Train########
train.model$varSelRF[['shrunken']] <-build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                              net.fea$varSelRF[c(3,4,5,6)],
                                                              stages.levels.comb)
train.model$varSelRF[['svm']] <-build.svm.classifier(vst_tumor_tum, tr.ind, 
                                                               net.fea$varSelRF[c(3,4,5,6)],
                                                               stages.levels.comb)
train.model$varSelRF[['rf']] <-build.rf.classifier(vst_tumor_tum, tr.ind, 
                                                               net.fea$varSelRF[c(3,4,5,6)],
                                                               stages.levels.comb)
train.model$varSelRF[['nb']] <-build.nb.classifier(vst_tumor_tum, tr.ind, 
                                                               net.fea$varSelRF[c(3,4,5,6)],
                                                               stages.levels.comb)
########CV########
cv.model$varSelRF[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, net.fea$varSelRF[c(3,4,5,6)],
                                               train.model$varSelRF$shrunken, tr.ind, 
                                               stages.levels.comb)
cv.model$varSelRF[['rf']] <- cv.rf.list(vst_tumor_tum, 10, net.fea$varSelRF[c(3,4,5,6)],
                                               tr.ind, stages.levels.comb)
cv.model$varSelRF[['svm']] <- cv.svm.list(vst_tumor_tum, 10, net.fea$varSelRF[c(3,4,5,6)],
                                        tr.ind, stages.levels.comb)
cv.model$varSelRF[['nb']] <- cv.nb.list(vst_tumor_tum, 10, net.fea$varSelRF[c(3,4,5,6)],
                                        tr.ind, stages.levels.comb)
cv.model$varSelRF[['knn']] <- cv.knn.list(vst_tumor_tum, 10, net.fea$varSelRF[c(3,4,5,6)],
                                        tr.ind, stages.levels.comb)
########Test#######
test.pred$varSelRF[['shrunken']] <- predict.shrunken(train.model$varSelRF$shrunken, 
                                                     net.fea$varSelRF[c(3,4,5,6)],
                                                     vst_tumor_tum[te.ind, ], 
                                                     cv.model$varSelRF$shrunken$thr)
test.pred$varSelRF[['rf']] <- predict.model(train.model$varSelRF$rf, 
                                            net.fea$varSelRF[c(3,4,5,6)],
                                            vst_tumor_tum[te.ind, ])
test.pred$varSelRF[['svm']] <- predict.model(train.model$varSelRF$svm, 
                                            net.fea$varSelRF[c(3,4,5,6)],
                                            vst_tumor_tum[te.ind, ])
test.pred$varSelRF[['nb']] <- predict.model(train.model$varSelRF$nb, 
                                             net.fea$varSelRF[c(3,4,5,6)],
                                             vst_tumor_tum[te.ind, ])
test.pred$varSelRF[['knn']] <- predict.knn(cv.model$varSelRF$knn, 
                                             net.fea$varSelRF[c(3,4,5,6)],
                                             vst_tumor_tum[tr.ind,],
                                             vst_tumor_tum[te.ind, ], 
                                           stages.levels.comb[tr.ind])
#############VarSelRF features###########

##############DeSeq2 features###########
train.model[['deseq2_2']] <- list()
cv.model[['deseq2_2']] <- list()
test.pred[['deseq2_2']] <- list()

##########Train###############
deseq2_2 <- lapply(net.fea$deseq2[c(3,4,5,6)], function(genes.list)
  {
  genes.list[['2 fold']]  
})
deseq2_1.5 <- lapply(net.fea$deseq2[c(3,4,5,6)], function(genes.list)
{
  genes.list[['1.5 fold']]  
})
deseq2_1 <- lapply(net.fea$deseq2[c(3,4,5,6)], function(genes.list)
{
  genes.list[['1 fold']]  
})
deseq2_2.5 <- lapply(net.features.updated$deseq2[c(3,4,5,6)], function(genes.list)
{
  genes.list[['2.5 fold']]  
})

deseq2_3 <- lapply(net.features.updated$deseq2[c(3,4,5,6)], function(genes.list)
{
  genes.list[['3 fold']]  
})



train.model$deseq2_2[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                               deseq2_2, stages.levels.comb)
train.model$deseq2_2[['svm']] <- build.svm.classifier(vst_tumor_tum, tr.ind, 
                                                                deseq2_2, stages.levels.comb)
train.model$deseq2_2[['nb']] <- build.nb.classifier(vst_tumor_tum, tr.ind, 
                                                                deseq2_2, stages.levels.comb)
train.model$deseq2_2[['rf']] <- build.rf.classifier(vst_tumor_tum, tr.ind, 
                                                                deseq2_2, stages.levels.comb)
train.model$deseq2_1.5[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                                deseq2_1.5, stages.levels.comb)
train.model$deseq2_1.5[['svm']] <- build.svm.classifier(vst_tumor_tum, tr.ind, 
                                                      deseq2_1.5, stages.levels.comb)
train.model$deseq2_1.5[['nb']] <- build.nb.classifier(vst_tumor_tum, tr.ind, 
                                                    deseq2_1.5, stages.levels.comb)
train.model$deseq2_1.5[['rf']] <- build.rf.classifier(vst_tumor_tum, tr.ind, 
                                                    deseq2_1.5, stages.levels.comb)
train.model$deseq2_1[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                                deseq2_1, stages.levels.comb)
train.model$deseq2_1[['svm']] <- build.svm.classifier(vst_tumor_tum, tr.ind, 
                                                      deseq2_1, stages.levels.comb)
train.model$deseq2_1[['nb']] <- build.nb.classifier(vst_tumor_tum, tr.ind, 
                                                    deseq2_1, stages.levels.comb)
train.model$deseq2_1[['rf']] <- build.rf.classifier(vst_tumor_tum, tr.ind, 
                                                    deseq2_1, stages.levels.comb)

##########CV###################
cv.model$deseq2_2[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_2,
                                              train.model$deseq2_2$shrunken, tr.ind, 
                                               stages.levels.comb)
cv.model$deseq2_1.5[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_1.5,
                                               train.model$deseq2_1.5$shrunken, tr.ind, 
                                               stages.levels.comb)
cv.model$deseq2_1[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_1,
                                               train.model$deseq2_1$shrunken, tr.ind, 
                                               stages.levels.comb)
cv.model$deseq2_2.5[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_2.5,
                                          train.model$deseq2_2.5$shrunken, train.indexes, 
                                               stages.levels.comb)


cv.model$deseq2_2[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_2,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_2[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_2,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_2[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_2,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_2[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_2,
                                          tr.ind, stages.levels.comb)

cv.model$deseq2_1.5[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_1.5,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_1.5[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_1.5,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_1.5[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_1.5,
                                            tr.ind, stages.levels.comb)
cv.model$deseq2_1.5[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_1.5,
                                            tr.ind, stages.levels.comb)

cv.model$deseq2_1[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_1,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_1[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_1,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_1[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_1,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_1[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_1,
                                          tr.ind, stages.levels.comb)

cv.model$deseq2_2.5[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_2.5,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_2.5[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_2.5,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_2.5[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_2.5,
                                            tr.ind, stages.levels.comb)
cv.model$deseq2_2.5[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_2.5,
                                            tr.ind, stages.levels.comb)


cv.model$deseq2_3[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_3,
                                          train.indexes, stages.levels.comb)
cv.model$deseq2_3[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_3,
                                          train.indexes, stages.levels.comb)
cv.model$deseq2_3[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_3,
                                            train.indexes, stages.levels.comb)
cv.model$deseq2_3[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_3,
                                            train.indexes, stages.levels.comb)

##########Test###############
test.pred$deseq2_2[['shrunken']] <- predict.shrunken(train.model$deseq2_2$shrunken, 
                                                     deseq2_2,
                                                     vst_tumor_tum[te.ind, ], 
                                                     cv.model$deseq2_2$shrunken$thr)
test.pred$deseq2_2[['svm']] <- predict.model(train.model$deseq2_2$svm, 
                                                     deseq2_2,
                                                     vst_tumor_tum[te.ind, ])
test.pred$deseq2_2[['rf']] <- predict.model(train.model$deseq2_2$rf, 
                                             deseq2_2,
                                             vst_tumor_tum[te.ind, ])
test.pred$deseq2_2[['nb']] <- predict.model(train.model$deseq2_2$nb, 
                                             deseq2_2,
                                             vst_tumor_tum[te.ind, ])
test.pred$deseq2_2[['knn']] <- predict.knn(cv.model$deseq2_2$knn, 
                                    deseq2_2, vst_tumor_tum[tr.ind,],
                                    vst_tumor_tum[te.ind, ], stages.levels.comb[tr.ind])

test.pred$deseq2_1.5[['shrunken']] <- predict.shrunken(train.model$deseq2_1.5$shrunken, 
                                                     deseq2_1.5,
                                                     vst_tumor_tum[te.ind, ], 
                                                     cv.model$deseq2_1.5$shrunken$thr)
test.pred$deseq2_1.5[['svm']] <- predict.model(train.model$deseq2_1.5$svm, 
                                             deseq2_1.5,
                                             vst_tumor_tum[te.ind, ])
test.pred$deseq2_1.5[['rf']] <- predict.model(train.model$deseq2_1.5$rf, 
                                            deseq2_1.5,
                                            vst_tumor_tum[te.ind, ])
test.pred$deseq2_1.5[['nb']] <- predict.model(train.model$deseq2_1.5$nb, 
                                            deseq2_1.5,
                                            vst_tumor_tum[te.ind, ])
test.pred$deseq2_1.5[['knn']] <- predict.knn(cv.model$deseq2_1.5$knn, 
                                           deseq2_1.5, vst_tumor_tum[tr.ind,],
                                           vst_tumor_tum[te.ind, ], stages.levels.comb[tr.ind])

test.pred$deseq2_1[['shrunken']] <- predict.shrunken(train.model$deseq2_1$shrunken, 
                                                     deseq2_1,
                                                     vst_tumor_tum[te.ind, ], 
                                                     cv.model$deseq2_1$shrunken$thr)
test.pred$deseq2_1[['svm']] <- predict.model(train.model$deseq2_1$svm, 
                                             deseq2_1,
                                             vst_tumor_tum[te.ind, ])
test.pred$deseq2_1[['rf']] <- predict.model(train.model$deseq2_1$rf, 
                                            deseq2_1,
                                            vst_tumor_tum[te.ind, ])
test.pred$deseq2_1[['nb']] <- predict.model(train.model$deseq2_1$nb, 
                                            deseq2_1,
                                            vst_tumor_tum[te.ind, ])
test.pred$deseq2_1[['knn']] <- predict.knn(cv.model$deseq2_1$knn, 
                                           deseq2_1, vst_tumor_tum[tr.ind,],
                                           vst_tumor_tum[te.ind, ], stages.levels.comb[tr.ind])

#####SAM Features###########
sam.fea <- list()
sam.fea[['2 fold']] <- lapply(net.fea$sam[c(3,4,5,6)], function(genes.list)
{
  genes.list[['2 fold']]  
})
sam.fea[['1.5 fold']] <- lapply(net.fea$sam[c(3,4,5,6)], function(genes.list)
{
  genes.list[['1.5 fold']]  
})
sam.fea[['1 fold']]  <- lapply(net.fea$sam[c(3,4,5,6)], function(genes.list)
{
  genes.list[['1 fold']]  
})

train.model$sam_2[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                                sam.fea$`2 fold`, stages.levels.comb)
train.model$sam_2[['svm']] <- build.svm.classifier(vst_tumor_tum, tr.ind, 
                                                   sam.fea$`2 fold`, stages.levels.comb)
train.model$sam_2[['nb']] <- build.nb.classifier(vst_tumor_tum, tr.ind, 
                                                 sam.fea$`2 fold`, stages.levels.comb)
train.model$sam_2[['rf']] <- build.rf.classifier(vst_tumor_tum, tr.ind, 
                                                 sam.fea$`2 fold`, stages.levels.comb)
train.model$sam_1.5[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                               sam.fea$`1.5 fold`, stages.levels.comb)
train.model$sam_1.5[['svm']] <- build.svm.classifier(vst_tumor_tum, tr.ind, 
                                                     sam.fea$`1.5 fold`, stages.levels.comb)
train.model$sam_1.5[['nb']] <- build.nb.classifier(vst_tumor_tum, tr.ind, 
                                                   sam.fea$`1.5 fold`, stages.levels.comb)
train.model$sam_1.5[['rf']] <- build.rf.classifier(vst_tumor_tum, tr.ind, 
                                                   sam.fea$`1.5 fold`, stages.levels.comb)
train.model$sam_1[['shrunken']] <- build.shrunken.classifier(vst_tumor_tum, tr.ind, 
                                                             sam.fea$`1 fold`, stages.levels.comb)
train.model$sam_1[['svm']] <- build.svm.classifier(vst_tumor_tum, tr.ind, 
                                                   sam.fea$`1 fold`, stages.levels.comb)
train.model$sam_1[['nb']] <- build.nb.classifier(vst_tumor_tum, tr.ind, 
                                                 sam.fea$`1 fold`, stages.levels.comb)
train.model$sam_1[['rf']] <- build.rf.classifier(vst_tumor_tum, tr.ind, 
                                                    sam.fea$`1 fold`, stages.levels.comb)

##########CV###################
cv.model$sam_2[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_2,
                                               train.model$deseq2_2$shrunken, tr.ind, 
                                               stages.levels.comb)
cv.model$sam_1.5[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_1.5,
                                                 train.model$deseq2_1.5$shrunken, tr.ind, 
                                                 stages.levels.comb)
cv.model$sam_1[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, deseq2_1,
                                               train.model$deseq2_1$shrunken, tr.ind, 
                                               stages.levels.comb)

cv.model$deseq2_2[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_2,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_2[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_2,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_2[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_2,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_2[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_2,
                                          tr.ind, stages.levels.comb)

cv.model$deseq2_1.5[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_1.5,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_1.5[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_1.5,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_1.5[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_1.5,
                                            tr.ind, stages.levels.comb)
cv.model$deseq2_1.5[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_1.5,
                                            tr.ind, stages.levels.comb)

cv.model$deseq2_1[['rf']] <- cv.rf.list(vst_tumor_tum, 10, deseq2_1,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_1[['nb']] <- cv.nb.list(vst_tumor_tum, 10, deseq2_1,
                                        tr.ind, stages.levels.comb)
cv.model$deseq2_1[['svm']] <- cv.svm.list(vst_tumor_tum, 10, deseq2_1,
                                          tr.ind, stages.levels.comb)
cv.model$deseq2_1[['knn']] <- cv.knn.list(vst_tumor_tum, 10, deseq2_1,
                                          tr.ind, stages.levels.comb)


save(cv.model, file = 'environment/accuracy_feature/updated/new_data/cv_model.RData')
save(test.pred, file = 'environment/accuracy_feature/updated/new_data/test_pred.RData')
save(train.model, file = 'environment/accuracy_feature/updated/new_data/train_model.RData')
