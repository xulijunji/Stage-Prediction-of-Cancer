install.packages('~/varSelRF_0.7-5.tar.gz', repos = NULL, type = 'source')
library(varSelRF)
library(snow)


df.stages$sample.id == colnames(exp_fpqm_tumor)

stages.levels <- df.stages$stage[indexes.stages.reported]
stages.levels <- droplevels(stages.levels)


cl <- makeMPIcluster(8)
tumor.fpqm.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels, 
                                 keep.forest = T)
tumor.fpqm.boot <- varSelRFBoot(t(exp_fpqm_tumor_reported), stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.fpqm.varSelRF)
tumor.fpqm.log.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_log_reported), Class = stages.levels, 
                                keep.forest = T)
tumor.fpqm.log.boot <- varSelRFBoot(t(exp_fpqm_tumor_log_reported), stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.fpqm.log.varSelRF)
tumor.exp.varSelRF <- varSelRF(xdata = t(exp_prof_tumor_reported), Class = stages.levels, 
                                keep.forest = T)
tumor.exp.boot <- varSelRFBoot(t(exp_prof_tumor_reported), stages.levels, bootnumber = 10, usingCluster = T,
                                TheCluster = cl, srf= tumor.exp.varSelRF)

tumor.fpqm.comb.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_reported), Class = stages.levels.comb, 
                                     keep.forest = T)
tumor.fpqm.comb.boot <- varSelRFBoot(t(exp_fpqm_tumor_reported), stages.levels.comb, bootnumber = 10, usingCluster = T,
                                     TheCluster = cl, srf= tumor.fpqm.comb.varSelRF)
tumor.fpqm.log.comb.varSelRF <- varSelRF(xdata = t(exp_fpqm_log_tumor_reported), Class = stages.levels.comb, 
                                         keep.forest = T)
tumor.fpqm.log.comb.boot <- varSelRFBoot(t(exp_fpqm_log_tumor_reported), stages.levels.comb, bootnumber = 10, usingCluster = T,
                                         TheCluster = cl, srf= tumor.fpqm.log.comb.varSelRF)

tumor.fpqm.over.varSelRF <- varSelRF(xdata = t(exp_fpqm_tumor_reported_over), Class = stages.oversamp,
                                     keep.forest = T)
tumor.fpqm.over.boot <- varSelRFBoot(t(exp_fpqm_tumor_reported_over), Class = stages.oversamp, bootnumber = 10, usingCluster = T,
                                     TheCluster = cl, srf= tumor.fpqm.over.varSelRF)

tumor.fpqm.smt.comb.varSelRF <-  varSelRF(xdata = data.frame(smt.comb[, -which(c('comb') == colnames(smt.comb))]),
                                          Class = smt.comb$comb, keep.forest = T)
tumor.nt.smt.comb.varSelRF <- varSelRF(xdata = data.frame(smt.comb.nt[, -which(c('comb') == colnames(smt.comb.nt))]),
                                            Class = smt.comb.nt$comb, keep.forest = T)
tumor.fpqm.under.comb.varSelRF <- varSelRF(xdata = und.comb.fpqm$X, Class = und.comb.fpqm$Y, keep.forest = T)
tumor.nt.under.comb.varSelRF <- varSelRF(xdata = und.comb.nt$X, Class = und.comb.nt$Y, keep.forest = T)
 

save(tumor.fpqm.smt.comb.varSelRF, file = 'environment/tumor_fpqm_smote_comb_var.RData')
save(tumor.nt.smt.comb.varSelRF, file = 'environment/tumor_nt_smote_comb_var.RData')
save(tumor.fpqm.under.comb.varSelRF, file = 'environment/tumor_fpqm_under_comb_var.RData')
save(tumor.nt.under.comb.varSelRF, file = 'environment/tumor_nt_under_comb_var.RData')
stopCluster(cl)


###
plot.lists <- function(varSelRF.lists, index)
{
  library(varSelRF)
  library(ggplot2)
  len <- length(varSelRF.lists)
  #print(len)
  par(mfrow=c(1,len))
  for(i in varSelRF.lists)
    plot(i,which=index)

}

plot.lists(list(tumor.fpqm.varSelRF, tumor.fpqm.log.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.varSelRF, tumor.fpqm.log.varSelRF), index = 2)

plot.lists(list(tumor.fpqm.comb.varSelRF, tumor.fpqm.log.comb.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.comb.varSelRF, tumor.fpqm.log.comb.varSelRF), index = 2)

plot.lists(list(tumor.fpqm.over.varSelRF), index = 1)
plot.lists(list(tumor.fpqm.over.varSelRF), index = 2)

plot(tumor.exp.boot, subject.names = df.stages$short.id)
selProbPlot(tumor.exp.boot)


dir.copy = paste(getwd(), 'results/tumor/varSelRF', sep = '/')

tumor.exp.varSelRF$selected.vars
tumor.fpqm.varSelRF$selected.vars
tumor.fpqm.log.varSelRF$selected.vars

tumor.exp.varSelRF$initialOrderedImportances[1:30]
tumor.fpqm.varSelRF$initialOrderedImportances[1:30]
tumor.fpqm.log.varSelRF$initialOrderedImportances[1:30]

##BootStrap
tumor.exp.boot$bootstrap.pred.error
tumor.fpqm.boot$bootstrap.pred.error
tumor.fpqm.log.boot$bootstrap.pred.error



tumor.exp.boot$all.data.randomForest$importance
tumor.fpqm.boot$all.data.randomForest$importance
tumor.fpqm.log.boot$all.data.randomForest$importance

View(tumor.exp.boot$class.predictions)

tumor.exp.boot$number.of.vars
tumor.fpqm.boot$number.of.vars
tumor.fpqm.log.boot$number.of.vars

tumor.exp.boot$overlap
tumor.fpqm.boot$overlap
tumor.fpqm.log.boot$overlap

save(tumor.exp.boot, file = 'environment/varSelRF_exp_boot.RData')
save(tumor.exp.varSelRF, file = 'environment/varSelRF_exp.RData')

save(tumor.fpqm.log.boot, file = 'environment/varSelRF_fpqm_log_boot.RData')
save(tumor.fpqm.log.varSelRF, file = 'environment/varSelRF_fpqk_log.RData')

save(tumor.fpqm.boot, file = 'environment/varSelRF_fpqm_boot.RData')
save(tumor.fpqm.varSelRF, file = 'environment/varSelRF_fpqm.RData')

save(tumor.fpqm.log.comb.boot, file = 'environment/tumor_fpqm_log_comb_boot.RData')
save(tumor.fpqm.log.comb.varSelRF, file = 'environment/tumor_fpqm_log_comb_varSelRF.RData')
save(tumor.fpqm.comb.boot, file = 'environment/tumor_fpqm_comb_boot.RData')
save(tumor.fpqm.comb.varSelRF, file = 'environment/tumor_fpqm_comb_varSelRF.RData')


save(tumor.indexes, file = 'environment/tumor_indexes.RData')
save(diff.ids, file = 'environment/diff_ids.RData') ##Used for tumor-control
save(indexes.stages.reported, file = 'environment/indexes_stages_reported.RData') ##Contain indexes for 
#which stage has been reported once the tumor.indexes are known

save(exp_prof, file = 'environment/exp_prof.RData' )
save(exp_fpqm, file = 'environment/exp_fpqm.RData' )


###Combining Stages
stages.levels.comb <- sapply(as.character(stages.levels), function(x)
  {
  if(x == 'stage i' || x == 'stage ii')
    x = 'stage i'
  else
    x = 'stage iv'
})
save(stages.levels.comb, file = 'environment/stages.level.comb.RData' )
