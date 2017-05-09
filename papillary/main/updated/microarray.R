source('main/updated/initialisation.R')
load('environment/accuracy_feature/updated/microarray_df.RData')
load('environment/accuracy_feature/updated/col_micro.RData')

col.indexes <- order(colnames(merged_final)[2:35])
micr_df <- merged_final[,c(2:35,38)]
micr_df <- micr_df[,c(col.indexes, 35)]
rownames(micr_df) <- micr_df$ens_id
micr_df <- micr_df[,-35]
micr_df <- t(micr_df)
micr_genes <- colnames(micr_df)

g_micr <- intersect(net.features.updated$shrunken$atleast_1, micr_genes)
g_micr <- intersect(g_micr, micr_genes)
res_micr <- final.res(micr_df, train.ind = c(1:34), test.ind = c(1:34), as.factor(sample_micro_info$class), 
                      g_micr, 10)

res.micr.class <- lapply(list(net.features.updated$shrunken[c(3,4,5,6)],
                              net.features.updated$varSelRF[c(3,4,5,6)],
                              deseq2_1, deseq2_1.5, deseq2_2), function(genes.list)
                                {
                                change.microarray(micr_df, genes.list, micr_genes, 
                                                  sample_micro_info$class, 10)
                              })
names(res.micr.class) <- c('shrunken', 'varSelRF', 'deseq2_1fold', 'deseq2_1.5fold', 'deseq2_2fold')

res.micr.stage <- lapply(list(net.features.updated$shrunken[c(3,4,5,6)],
                              net.features.updated$varSelRF[c(3,4,5,6)],
                              deseq2_1, deseq2_1.5, deseq2_2), function(genes.list)
                              {
                                change.microarray(micr_df, genes.list, micr_genes, 
                                                  sample_micro_info$stage, 10)
                              })
names(res.micr.stage) <- c('shrunken', 'varSelRF', 'deseq2_1fold', 'deseq2_1.5fold', 'deseq2_2fold')

micr.dfs.class <- lapply(res.micr.class, function(x)
{
  publish.results(x)
})
micr.dfs.stage <- lapply(res.micr.stage, function(x)
{
  publish.results(x)
})
micr.dfs.stage$varSelRF <- publish.results(res.micr.stage$varSelRF, indexes = c(1,2,3))
get.conf.mat(res_micr[[1]])
get.aucs(res$nb)
get.aucs(res_micr[[1]])
