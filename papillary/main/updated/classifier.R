source('main/updated/initialisation.R')

train.indexes <- sort(unlist(gr.train))
train.model <- list()
cv.model <- list()

##############Shrunken features###############
train.model[['shrunken']] <- list()
cv.model[['shrunken']] <- list()
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

######Shrunken CV######
cv.model$shrunken[['shrunken']] <- cv.shrunken(vst_tumor_tum, 10, net.features.updated$shrunken[c(3,4,5,6)],
                                               train.model$shrunken$shrunken, train.indexes, 
                                               stages.levels.comb)
length(train.model$shrunken$shrunken)
