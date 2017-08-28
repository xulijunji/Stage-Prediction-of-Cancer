setwd('Honours/Stage-Prediction-of-Cancer/papillary/')
library(Biocomb)
library(class)
library(randomForest)
library(e1071)
library(DESeq2)
source('stability_gene_selection/helper_stable.R')
source('../function.R')
source('feature_sel/class_loop.R')
load('environment/stages_levels.RData')
load('environment/diff_genes.RData')
load('environment/stage.index.RData')
load('environment/stages.level.comb.RData')
load('environment/dds_tumor_reported_normal_stage.RData')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/vs_normal_tumor_repored.RData')
load('environment/filter_based_genes.RData_based.RData')
load('environment/filter_based_1_2.RData')

genes.most.varying <- rownames(vs.normal.tumor.reported)[get.imp.genes(3,
                                                                       assay(vs.normal.tumor.reported), 10000)]
sum(rownames(sample.info.all.rep) == colnames(vs.normal.tumor.reported)) == length(rownames(sample.info.all.rep))
indexes.to.remove <- which(sample.info.all.rep$type == 'N')
vs.normal.tumor.reported.varying <- t(assay(vs.normal.tumor.reported)[genes.most.varying,
                                                                      -indexes.to.remove])
sum(rownames(vs.normal.tumor.reported.varying) == sample.info.all.rep$sample.names[sample.info.all.rep$type == 'T'])

rf.relief <- randomForest(vs.normal.tumor.reported.varying[,
                                              as.character(filter.based.genes$relief$Biomarker[1:300])],
                          y = stages.levels)
rf.relief$confusion
stage.ind.1.2 <- sort(unlist(stage.ind[c(1,2)]))
stage.ind.1.3 <- sort(unlist(stage.ind[c(1,3)]))

filter.based.1.2 <- list()
filter.based.1.2[['FastFilter']] <- select.fast.filter(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                                                             as.character(stages.levels[stage.ind.1.2])), 
                                                       disc.method = 'MDL', threshold = 0.05)
filter.based.1.2[['InformationGain']] <- select.inf.gain(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                                                                  as.character(stages.levels[stage.ind.1.2])), 
                                                            disc.method = 'MDL')
filter.based.1.2[['SymmUnc']] <- select.inf.symm(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                                                                  as.character(stages.levels[stage.ind.1.2])), 
                                                            disc.method = 'MDL')
filter.based.1.2[['Relief']] <- select.relief(as.data.frame(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                                                    as.character(stages.levels[stage.ind.1.2]))))

filter.based.1.3 <- list()
filter.based.1.2[['FastFilter']] <- select.fast.filter(cbind(vs.normal.tumor.reported.varying[stage.ind.1.3,],
                                                             as.character(stages.levels[stage.ind.1.3])), 
                                                       disc.method = 'MDL', threshold = 0.05)

rf.1.2 <- list()
rf.1.2[['filter']] <- list()
rf.1.2$filter[['symmUnc']] <- randomForest(vs.normal.tumor.reported.varying[stage.ind.1.2,
                      filter.based.1.2$SymmUnc$Biomarker[1:300]], droplevels(stages.levels[stage.ind.1.2]))


rf.1.3 <- list()
rf.1.3[['filter']] <- list()
rf.1.3[['filter']][['fast_filt']] <- randomForest(vs.normal.tumor.reported.varying[stage.ind.1.3,
                                                       filter.based.1.2$SymmUnc$Biomarker[1:300]],
                                                  droplevels(stages.levels[stage.ind.1.3]))







class.loop.1.2 <- list()
class.loop.1.2[['InformationGain']] <- select.fast.filter(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                                                                as.character(stages.levels[stage.ind.1.2])), 
                                                          disc.method = 'MDL', threshold = 0.05)
  classifier.loop(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                        as.character(stages.levels[stage.ind.1.2])), 
                  feature.selection=c("InformationGain"), threshold = 0.1, no.feat = 300,
                  method.cross = 'fold-crossval')


class.loop.1.2[["symmetrical.uncertainty"]] <-
  classifier.loop(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                        as.character(stages.levels[stage.ind.1.2])), 
                  feature.selection=c("symmetrical.uncertainty"), threshold = 0.1, no.feat = 300)


class.loop.1.2[["FastFilter"]] <-
  classifier.loop(cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],
                        as.character(stages.levels[stage.ind.1.2])), 
                  feature.selection=c("FastFilter"), threshold = 0.1, no.feat = 300)

a <- cbind(vs.normal.tumor.reported.varying[stage.ind.1.2,],droplevels(stages.levels[stage.ind.1.2]))
colnames(a)[10001] = 'class.label'
svm.mod <- svm(class.label~.,data = a[,c(1:10, 10001)])            
