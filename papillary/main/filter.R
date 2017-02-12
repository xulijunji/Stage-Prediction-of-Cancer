library(DESeq2)
library(Biocomb)
source('decomp.R')
source('main/final.R')

load('environment/dds_nor_tum_comb.RData')
load('environment/first_trial_shrunken_classifier.RData')
vs_normal_comb_reported <- t(assay(vst(dds_nor_tum_comb)))
genes.most.varying <- colnames(vs_normal_comb_reported)[get.imp.genes(3,
                                                                      t(vs_normal_comb_reported), 10000)]
tumor.ind.vs <- which(colData(dds_nor_tum_comb)[,2] == 'T')
match.ind <- match(colData(dds_nor_tum_comb)[,1][tumor.ind.vs],
                   df.stage.tumor.rep$sample.id)

sum(droplevels(colData(dds_nor_tum_comb)[tumor.ind.vs,4]) ==
      stages.levels.comb)
vs_normal_comb_reported_var <- vs_normal_comb_reported[tumor.ind.vs, genes.most.varying]

gr <- first.trial$gr

features.filter.based <- list()
features.filter.based[['ffs']] <- do.fast.filter(gr, vs_normal_comb_reported_var, stages.levels.comb)
features.filter.based[['forCorr']] <- do.forward.corr(gr, vs_normal_comb_reported_var, stages.levels.comb)

View(features.filter.based$ffs[[1]])
save(features.filter.based, file = 'environment/feautures_filter.RData')
load('environment/feautures_filter.RData')


