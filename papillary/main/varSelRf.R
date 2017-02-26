###Using varselRf for feature selection
library(varSelRF)
library(DESeq2)

source('decomp.R')
source('main/final.R')
source('varselRF/varSelRf_altered.R')
load('environment/first_trial_shrunken_classifier.RData')
load('environment/feature_varselRF.RData')
load('environment/features_varSelRF1_comb.RData')
load('environment/accuracy_feature/net_features.RData')
load('environment/accuracy_feature/classifer_list.RData')
load('environment/dds_nor_tum_comb.RData')
##Using vst_normal_reported,tumor.ind.vs from main.R(shrunken)


gr <- first.trial$gr
features.varSelRF <- do.varselRF(gr, vst_normal_reported[tumor.ind.vs,],
                                 stages.levels.comb)
features.varSelRF1 <- do.varselRF(gr, vs_normal_comb_reported[tumor.ind.vs,],
                                 stages.levels.comb)
save(features.varSelRF1, file = 'environment/features_varSelRF1_comb.RData')

genes.varSelRf.list <- get.varSelRf.genes(features.varSelRF)
genes.varSelRf.list1 <- get.min.oob.varselRf(net.features$varSelRf$genes.object)

net.features[['varSelRf']] <- list()
net.features$varSelRf[['genes.object']] <- features.varSelRF1
net.features$varSelRf[['genes.list']] <- get.min.oob.varselRf(net.features$varSelRf$genes.object)

sapply(genes.varSelRf.list1, length)
length(intersect(genes.varSelRf.list1[[2]], genes.varSelRf.list1[[1]]))

net.features$varSelRf[['atleast_1']] <- get.genes.common(genes.varSelRf.list1, 1)
net.features$varSelRf[['atleast_3']] <- get.genes.common(genes.varSelRf.list1, 3)
net.features$varSelRf[['atleast_5']] <- get.genes.common(genes.varSelRf.list1, 5)
length(net.features$varSelRf$atleast_1)
length(intersect(net.features$varSelRf$atleast_5, net.features$shrunken$atleast_5))

classifier.list$varSelRf <- list()
classifier.list$varSelRf[['atleast_1']] <- do.rf(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$varSelRf$atleast_1,
                                                 stages.levels.comb, list(), 
                                                 list())
classifier.list$varSelRf[['atleast_1']] <- do.knn(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$varSelRf$atleast_1,
                                                 stages.levels.comb, classifier.list$varSelRf$atleast_1[[1]], 
                                                 classifier.list$varSelRf$atleast_1[[2]])
classifier.list$varSelRf[['atleast_1']] <- do.naive(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$varSelRf$atleast_1,
                                                 stages.levels.comb, classifier.list$varSelRf$atleast_1[[1]], 
                                                 classifier.list$varSelRf$atleast_1[[2]])
classifier.list$varSelRf[['atleast_1']] <- do.svm(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$varSelRf$atleast_1,
                                                 stages.levels.comb, classifier.list$varSelRf$atleast_1[[1]], 
                                                 classifier.list$varSelRf$atleast_1[[2]])

classifier.list$varSelRf[['atleast_3']] <- list()
classifier.list$varSelRf[['atleast_3']] <- do.rf(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$varSelRf$atleast_3,
                                                 stages.levels.comb, list(), 
                                                 list())
classifier.list$varSelRf[['atleast_3']] <- do.knn(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$varSelRf$atleast_3,
                                                 stages.levels.comb, classifier.list$varSelRf$atleast_3[[1]], 
                                                 classifier.list$varSelRf$atleast_3[[2]])
classifier.list$varSelRf[['atleast_3']] <- do.svm(first.trial$gr, 
                                                  vs_normal_comb_reported[tumor.ind.vs,], 
                                                  net.features$varSelRf$atleast_3,
                                                  stages.levels.comb, classifier.list$varSelRf$atleast_3[[1]], 
                                                  classifier.list$varSelRf$atleast_3[[2]])
classifier.list$varSelRf[['atleast_3']] <- do.naive(first.trial$gr, 
                                                  vs_normal_comb_reported[tumor.ind.vs,], 
                                                  net.features$varSelRf$atleast_3,
                                                  stages.levels.comb, classifier.list$varSelRf$atleast_3[[1]], 
                                                  classifier.list$varSelRf$atleast_3[[2]])

classifier.list$varSelRf[['atleast_5']] <- list()
classifier.list$varSelRf[['atleast_5']] <- do.knn(first.trial$gr, 
                                                  vs_normal_comb_reported[tumor.ind.vs,], 
                                                  net.features$varSelRf$atleast_5,
                                                  stages.levels.comb, list(), 
                                                  list())
classifier.list$varSelRf[['atleast_5']] <- do.rf(first.trial$gr, 
                                                  vs_normal_comb_reported[tumor.ind.vs,], 
                                                  net.features$varSelRf$atleast_5,
                                                  stages.levels.comb, classifier.list$varSelRf$atleast_5[[1]], 
                                                  classifier.list$varSelRf$atleast_5[[2]])
classifier.list$varSelRf[['atleast_5']] <- do.naive(first.trial$gr, 
                                                 vs_normal_comb_reported[tumor.ind.vs,], 
                                                 net.features$varSelRf$atleast_5,
                                                 stages.levels.comb, classifier.list$varSelRf$atleast_5[[1]], 
                                                 classifier.list$varSelRf$atleast_5[[2]])
classifier.list$varSelRf[['atleast_5']] <- do.svm(first.trial$gr, 
                                                    vs_normal_comb_reported[tumor.ind.vs,], 
                                                    net.features$varSelRf$atleast_5,
                                                    stages.levels.comb, classifier.list$varSelRf$atleast_5[[1]], 
                                                    classifier.list$varSelRf$atleast_5[[2]])
save(classifier.list, file = 'environment/accuracy_feature/classifer_list.RData')
save(net.features, file = 'environment/accuracy_feature/net_features.RData')
