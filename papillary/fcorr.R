library(Biocomb)
library(class)
library(randomForest)
library(e1071)
source('../stability_gene_selection/helper_stable.R')

load('environment/stages_levels.RData')
load('environment/diff_genes.RData')
load('environment/filter_based_genes.RData')
load('environment/only_tumor_reported.RData')

filter.based.genes = list()
data.fpqm.filter = t(only.tumor.reported$dfs$fpqm)
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

pred.knn.fpqm.fastcorr <- knn.cv(data.fpqm.filter[,filter.based.genes$ffilter$Biomarker[1:6]], cl = stages.levels, k = 3)
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

gamma = 0, kernel = 'linear', cost = 1,
class.weights =if(length(levels(stages.levels)) == 4) c('stage i' = 1, 'stage ii' =1, 
                                                        'stage iii' = 1, 'stage iv' =1)
else c('stage i' = 1, 'stage iv' = 1)