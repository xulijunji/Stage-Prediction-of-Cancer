library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query = query, directory = '~/ccRCC/')
data <- GDCprepare(query, directory = '~/ccRCC/')
save(data, file = 'environment/kirc_data.RData')

load('../environment/kirc_data.RData')
table(colData(data)[,4])
tumor.inds <- which(colData(data)[,4] == 'TAP' | colData(data)[,4] == 'TP')
ccrcc.tum.data <- data[, tumor.inds]
not.rep.inds <- which(colData(ccrcc.tum.data)[,10] == 'not reported')
ccrcc.tum.data <- ccrcc.tum.data[,-not.rep.inds]
table(colData(ccrcc.tum.data)[,10])
stages.ccrcc <- colData(ccrcc.tum.data)[,10]
stages.ccrcc.comb <- sapply(stages.ccrcc, function(stage)
{
  if(stage == 'stage i' | stage == 'stage ii')
    'stage i'
  else
    'stage iv'
})

source('main/updated/classifier.R')
load('environment/accuracy_feature/updated/new_data/train_model.RData')
load('environment/accuracy_feature/updated/new_data/fea.RData')
load('environment/accuracy_feature/updated/new_data/cv_model.RData')
source('decomp.R')
source('main/updated/helper_func.R')
source('main/final.R')
library(DESeq2)
norm.data <- t(vst(assay(ccrcc.tum.data)))
fea.trial.list.crc <- lapply(fea.trial.list, function(x) 
{lapply(x, function(x)intersect(x, colnames(norm.data)))})

load('environment/accuracy_feature/updated/vst_tum_rep.RData')

prcc.common.genes <- intersect(colnames(vst_tumor_tum), colnames(norm.data))
comb.data <- rbind(vst_tumor_tum[,prcc.common.genes], norm.data[, prcc.common.genes])
pr <- get.test.pred(data = t(norm.data), te.ind = seq(ncol(norm.data)), 
                    fea.list = fea.trial.list.crc[c(1,2)], stages = NULL, tr.ind = 1:260,
                    tr.model = train.trial.model, cv.model = cv.trial.model
                      )

ccrcc.pred <- get.test.pred(data = comb.data, te.ind = c(261:796), 
                    fea.list = fea.trial.list.crc[c(1,2,3,4,5,7,8)], stages = stages.levels.comb, 
                    tr.ind = 1:260,
                    tr.model = train.trial.model, cv.model = cv.trial.model
)

ccrcc.results <- get.test.res(ccrcc.pred, seq(length(stages.ccrcc.comb)), stages.ccrcc.comb, names(ccrcc.pred))

