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

g_micr <- intersect(net.features.updated$varSelRF$atleast_4, micr_genes)
g_micr <- intersect(g_micr, micr_genes)



res_micr <- final.res(micr_df, train.ind = c(1:34), test.ind = c(1:34), as.factor(sample_micro_info$stage), 
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

######Proper#####
load('environment/accuracy_feature/updated/new_data/fea.RData')

fea.micr <- list()
for(i in seq_along(fea.trial.list))
  fea.micr[[names(fea.trial.list)[i]]] <- lapply(fea.trial.list[[names(fea.trial.list)[i]]], function(g) intersect(g, micr_genes))

train.micr <- get.train.model(micr_df, tr.ind = c(1:34), fea.micr, as.factor(sample_micro_info$stage))
cv.micr <- get.cv.model(micr_df, 1:34, fea.micr, as.factor(sample_micro_info$stage), train.micr)
cv.micr.res <- get.cv.res(1:34, cv.micr, as.factor(sample_micro_info$stage), names(fea.micr))

micr.net.df <-create.net.df(cv.micr.res)
ggplot(micr.net.df, aes(x=Classifier, y=AUC, fill=Feature_Selection)) +
  geom_boxplot()

create.net.plot(micr.net.df)
get.aucs(cv.micr.res$`sam1 fold`$shrunken)

fea.old.list <- get.class.fea(net.features.updated)
fea.micr.old <- list()
for(i in seq_along(fea.old.list))
  fea.micr.old[[names(fea.old.list)[i]]] <- lapply(fea.old.list[[names(fea.old.list)[i]]], function(g) intersect(g, micr_genes))

micr.int <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage), 
                      intersect(micr_genes, g1), 10)
micr.wcgna <- final.res(micr_df, 1:34, 1:34, as.factor(sample_micro_info$stage), 
                        intersect(micr_genes, wcgna.genes.ent.ens), 10)
get.aucs(micr.int[[1]])
