source('main/updated/initialisation.R')
load('environment/accuracy_feature/updated/microarray_df.RData')
load('environment/accuracy_feature/updated/col_micro.RData')

col.indexes <- order(colnames(merged_final)[2:35])
micr_df <- merged_final[,c(2:35,38)]
micr_df <- micr_df[,c(col.indexes, 35)]
rownames(micr_df) <- micr_df$ens_id
micr_df <- micr_df[,-35]
micr_genes <- colnames(micr_df)

g_micr <- intersect(net.features.updated$deseq2$atleast_2$`1 fold`, micr_genes)
g_micr <- intersect(g_1.5_1, micr_genes)
res_micr <- final.res(micr_df, train.ind = c(1:34), test.ind = c(1:34), as.factor(sample_micro_info$stage), 
                      g_micr, 10)
get.conf.mat(res_micr[[1]])
