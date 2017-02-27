source('decomp.R')
source('main/updated/helper_func.R')
source('main/final.R')
#library(pamr)

load('environment/stages.level.comb.RData')
load('environment/accuracy_feature/vs_nor_comb.RData')
load('environment/accuracy_feature/tumor_ind_vs.RData')
load('environment/first_trial_shrunken_classifier.RData')
load('environment/accuracy_feature/gr.RData')

vst_tumor_tum <- vs_normal_comb_reported[tumor.ind.vs, ]

gr.updated <- build.groups(260,5)
stage.dist <- get.stage.distribution(gr.updated, stages.levels.comb)

test.group <- gr.updated[[1]]
train.group <- unlist(gr.updated[-1])

gr.train <- gr.updated[-1] ##in terms of 260

net.features.updated <- list()
######Shrunken####
net.features.updated[['shrunken']] <- list()
shrunken.features.ob <- get.shrunken.features(gr.train, vst_tumor_tum, stages.levels.comb)
net.features.updated$shrunken[['genes.object']] <- shrunken.features.ob$genes
net.features.updated$shrunken[['genes']] <- get.genes.shrunken(net.features.updated$shrunken$genes.object)
net.features.updated$shrunken[['atleast_1']] <- get.genes.common(net.features.updated$shrunken$genes, 1)
net.features.updated$shrunken[['atleast_2']] <- get.genes.common(net.features.updated$shrunken$genes, 2)
net.features.updated$shrunken[['atleast_4']] <- get.genes.common(net.features.updated$shrunken$genes, 4)
net.features.updated$shrunken[['genes_stage']] <- get.shrunken.group.stage(net.features.updated$shrunken$genes.object)
net.features.updated$shrunken[['atleast_dfs']] <- get.shrunken.common.stage(net.features.updated$shrunken$genes_stage, c(1,2,4))

#####