###Using varselRf for feature selection
library(varSelRF)
library(DESeq2)

source('decomp.R')
source('varselRF/varSelRf_altered.R')
load('environment/first_trial_shrunken_classifier.RData')
##Using vst_normal_reported,tumor.ind.vs from main.R(shrunken)


gr <- first.trial$gr
features.varSelRF <- do.varselRF(gr, vst_normal_reported[tumor.ind.vs,],
                                 stages.levels.comb)
features.varSelRF1 <-do.varselRF(gr, vst_normal_reported[tumor.ind.vs,],
                                 stages.levels.comb)
intersect(features.varSelRF[[3]]$selected.vars, features.varSelRF[[5]]$selected.vars)
