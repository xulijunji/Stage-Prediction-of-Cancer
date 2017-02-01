library(DESeq2)
load('environment/stages.level.comb.RData')
load('papillary/environment/dds_tumor_reported_normal_stage.RData')
load('h_short/papillary/environment/dds_tumor_reported.RData')
load('h_short/papillary/environment/first_trial_shrunken_classifier.RData')
source('papillary/across_tumors.R')
source('h_short/papillary/shrunken/final.R')

sample.df <- colData(dds_tumor_reported_normal_stage)
sample.df[,4] <- sapply(as.character(sample.df[,4]), function(stage)
  {
  if(stage == 'stage ii')
    stage = 'stage i'
  else if(stage == 'stage iii')
    stage = 'stage iv'
  else
    stage = stage
})
dds_nor_tum_comb <- 
  DESeqDataSetFromMatrix(counts(dds_tumor_reported_normal_stage),
                         colData = sample.df, design = ~stage.type)
dds_nor_tum_comb <- dds_nor_tum_comb[rowSums(assay(dds_nor_tum_comb)) > 2]
dds_nor_tum_comb <- DESeq(dds_nor_tum_comb)

sample.df <- colData(dds_tumor_reported)
sample.df[,3] <- sapply(as.character(sample.df[,3]), function(stage)
  {
  if(stage == 'stage ii')
    stage = 'stage i'
  else if(stage == 'stage iii')
    stage = 'stage iv'
  else
    stage = stage
})
gr = first.trial$gr
dds_obj <- create.Deseq2(gr, counts(dds_tumor_reported), 
                         colData = sample.df)
res.train <- do.Deseq2(first.trial$gr, dds_nor_tum_comb)
