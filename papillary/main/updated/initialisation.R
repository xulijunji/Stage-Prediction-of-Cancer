#library(pamr)

load('environment/stages.level.comb.RData')
#load('environment/accuracy_feature/dds_nor_tum_comb.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')
load('environment/accuracy_feature/tumor_ind_vs.RData')
#load('environment/first_trial_shrunken_classifier.RData')
load('environment/accuracy_feature/updated/gr.RData')
load('environment/accuracy_feature/updated/net_features_updated.RData')
load('environment/accuracy_feature/updated/gr_trial.RData')

source('decomp.R')
source('main/updated/helper_func.R')
source('main/final.R')

vst_tumor_tum <- vs_normal_comb_reported[tumor.ind.vs, ]
remove(vs_normal_comb_reported)
#gr.updated <- build.groups(260,5)
stage.dist <- get.stage.distribution(gr.updated, stages.levels.comb)
gr.updated.train <- gr.updated[-1]

#gr.trial <- build.groups(260,5)
#gr.trial.train <-
#get.stage.distribution(gr.trial, stages.levels.comb)

gr.trial.train <- gr.trial[-5]
test.indexes <- gr.updated[[1]]
train.indexes <- sort(unlist(gr.updated[-1]))

train.trial.ind <- sort(unlist(gr.trial[-5]))
test.trial.ind <- gr.trial[[5]]
#save(gr.trial, file = 'environment/accuracy_feature/updated/gr_trial.RData')
