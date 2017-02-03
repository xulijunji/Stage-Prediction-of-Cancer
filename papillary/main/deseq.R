library(DESeq2)
load('environment/stages.level.comb.RData')
load('environment/dds_tumor_reported_normal_stage.RData')
load('environment/dds_tumor_reported.RData')
load('environment/first_trial_shrunken_classifier.RData')
load('environment/dds_object.RData')
source('across_tumors.R')
source('main/final.R')

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
save(dds_nor_tum_comb, file = 'environment/dds_nor_tum_comb.RData')

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
res.train <- do.Deseq2(dds_obj)
deseq.genes.list <- sapply(res.train, function(x)
  {
  get.genes(x, 1.5, 0.01, 0.01)
})
sapply(deseq.genes.list, length)
deseq.genes.common <- get.intersect.genes(deseq.genes.list, seq(5))
length(intersect(deseq.genes.common, g2))
