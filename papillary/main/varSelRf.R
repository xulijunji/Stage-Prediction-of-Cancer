###Using varselRf for feature selection
library(varSelRF)
library(DESeq2)

source('decomp.R')
source('main/final.R')
source('varselRF/varSelRf_altered.R')
load('environment/first_trial_shrunken_classifier.RData')
load('environment/feature_varselRF.RData')
load('environment/features_varSel.RData')
##Using vst_normal_reported,tumor.ind.vs from main.R(shrunken)


gr <- first.trial$gr
features.varSelRF <- do.varselRF(gr, vst_normal_reported[tumor.ind.vs,],
                                 stages.levels.comb)
features.varSelRF1 <-do.varselRF(gr, vs_normal_comb_reported[tumor.ind.vs,],
                                 stages.levels.comb)
save(features.varSelRF1, file = 'environment/features_varSelRF1_comb.RData')

genes.varSelRf.list <- get.varSelRf.genes(features.varSelRF)
genes.varSelRf.list1 <- get.min.oob.varselRf(features.varSelRF1)
sapply(genes.varSelRf.list1, length)
length(intersect(genes.varSelRf.list1[[2]], genes.varSelRf.list1[[1]]))
length(intersect(genes.varSelRf, g2))
genes.varSelRf <- get.genes.common(genes.varSelRf.list1, 5)
