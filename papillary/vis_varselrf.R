load('environment/tumor_fpqm_comb_boot.RData')
load('environment/tumor_fpqm_comb_varSelRF.RData')
load('environment/tumor_fpqm_log_comb_boot.RData')
load('environment/tumor_fpqm_log_comb_varSelRF.RData')
load('environment/tumor_fpqm_over_boot.RData')
load('environment/tumor_fpqm_over_varSelRF.RData')

library(varSelRF)
intersect(tumor.fpqm.comb.varSelRF$selected.vars,tumor.fpqm.log.comb.varSelRF$selected.vars)

save.confusion.mat <- function(boot, dir)
tumor.exp.boot$all.data.randomForest$confusion
tumor.fpqm.boot$all.data.randomForest$confusion
tumor.fpqm.log.boot$all.data.randomForest$confusion
tumor.fpqm.log.comb.boot$all.data.randomForest$confusion
tumor.fpqm.comb.boot$all.data.randomForest$confusion
tumor.fpqm.over.boot$all.data.randomForest$confusion

tumor.fpqm.over.boot$all.data.vars


