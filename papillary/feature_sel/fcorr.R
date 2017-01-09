library(Biocomb)
library(class)
library(randomForest)
library(e1071)
library(DESeq2)
source('stability_gene_selection/helper_stable.R')
source('../function.R')

load('environment/stages_levels.RData')
load('environment/diff_genes.RData')
load('environment/stage.index.RData')
load('environment/stages.level.comb.RData')
load('environment/dds_tumor_reported_normal_stage.RData')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/vs_normal_tumor_repored.RData')
load('environment/filter_based_genes.RData_based.RData')

vs.normal.tumor.reported <-
  varianceStabilizingTransformation(dds_tumor_reported_normal_stage)
genes.most.varying <- rownames(vs.normal.tumor.reported)[get.imp.genes(3,
                                        assay(vs.normal.tumor.reported), 10000)]
sum(rownames(sample.info.all.rep) == colnames(vs.normal.tumor.reported)) == length(rownames(sample.info.all.rep))
indexes.to.remove <- which(sample.info.all.rep$type == 'N')
vs.normal.tumor.reported.varying <- t(assay(vs.normal.tumor.reported)[genes.most.varying,
                                                                  -indexes.to.remove])
sum(rownames(vs.normal.tumor.reported.varying) == sample.info.all.rep$sample.names[sample.info.all.rep$type == 'T'])


filter.based.genes = list()
filter.based.genes[['fast_filter']] <- select.fast.filter(cbind(vs.normal.tumor.reported.varying, stages.levels), 
                                                          disc.method = 'MDL', threshold = 0.05)
filter.based.genes[['forw.Corr']] <- select.forward.Corr(cbind(vs.normal.tumor.reported.varying,
                                           stages.levels), disc.method = 'MDL')
filter.based.genes[['relief']] <- select.relief(as.data.frame(cbind(vs.normal.tumor.reported.varying,
                                                            stages.levels)))
filter.based.genes[['forw_Wrap']] <- select.forward.wrapper(as.data.frame(cbind(vs.normal.tumor.reported.varying,
                                                                    stages.levels)))

data.fpqm.filter = req.dfs$fpqm
data.fpqm.filter.diff1 <- data.fpqm.filter[,diff.genes[[1]]]
data.fpqm.filter$class = stages.levels
typeof(as.data.frame(data.fpqm.filter))
filter.based.genes[['cfs']] = select.cfs(data.fpqm.filter.diff)
filter.based.genes[['ffilter']] = select.fast.filter(data.fpqm.filter, 'MDL', 0.1,numeric())
                                                                                                                                                  filter.based.genes[['ffilter']] = select.fast.filter(data.fpqm.filter.diff, 'MDL', 0.05,numeric())
filter.based.genes[['fcor']] = select.forward.Corr(data.fpqm.filter.diff1, 'MDL', numeric())
save(filter.based.genes, file = 'environment/filter_based_genes.RData')
load('environment/filter_based_genes.RData')
intersect(filter.based.genes$cfs$Biomarker, filter.based.genes$ffilter$Biomarker)

pred.knn.fpqm.fastcorr <- knn.cv(data.fpqm.filter[,filter.based.genes$ffilter$Biomarker[1:76]], cl = stages.levels, k = 3)
create.mat.error(table(stages.levels, pred.knn.fpqm.fastcorr))

pred.knn.fpqm.cfs <- knn.cv(data.fpqm.filter[,filter.based.genes$cfs$Biomarker], cl = stages.levels, k = 5)
create.mat.error(table(stages.levels, pred.knn.fpqm.cfs))

rf.fpqm.fastcorr <- randomForest(x = data.fpqm.filter[,filter.based.genes$ffilter$Biomarker], 
                                 y = stages.levels, strata = stages.levels, sampsize = c(20,15,10,15))
rf.fpqm.fastcorr$confusion

rf.fpqm.cfs <- randomForest(x = data.fpqm.filter[,filter.based.genes$cfs$Biomarker], 
                           y = stages.levels, strata = stages.levels, sampsize = c(20,15,10,15))
rf.fpqm.cfs$confusion

svm.fastcorr <- cv.svm.leave.one.out(data.fpqm.filter[,filter.based.genes$ffilter$Biomarker], stages.levels,
                                     cost = 1, class.weights = c('stage i' = 1, 'stage ii' = 5, 'stage iii' = 3, 'stage iv' = 7))

svm.cfs <- cv.svm.leave.one.out(data.fpqm.filter[,filter.based.genes$cfs$Biomarker], stages.levels)
svm.cfs
svm.fastcorr <- cv.svm.leave.one.out(data.fpqm.filter[,filter.based.genes$ffilter$Biomarker], stages.levels)
svm.fastcorr
