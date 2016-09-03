###Tumor samples
##Uses function.R for all processing and visualisation
##In it we create data frames/dds objects for only tumor samples.
##We also separate them into only those for which we stage data.
##Env variables saved df.stages, dds_tumor, dds_tumor_reported, only.tumor, only.tumor.reported
##Proj variables - df.stage.tumor.rep, df.stages, indexes.stages.reported, pca_tumor_rep

library(DESeq2)
load('environment/exp_prof.RData')
load('environment/exp_fpqm.RData')


pat.stages = rep('NA', length(sample.ids))
indexes.stages = match(pats[tumor.indexes], names(pat_stages))
pat.stages[tumor.indexes] = pat_stages[indexes.stages] #from extract_stage.R
names(pat.stages) = sample.ids

df.stages = data.frame(sample.id=sample.ids[tumor.indexes], patient.id = pats[tumor.indexes], stage =
                         unlist(pat_stages[indexes.stages]), short.id = c(1:289))
rownames(df.stages) = c(1:289)

indexes.stages.reported <- which(df.stages$stage != 'not reported') ##Contains indexes inside tumor.indexes

exp_prof_tumor = exp_prof[,tumor.indexes]
exp_prof_tumor_reported = exp_prof_tumor[,indexes.stages.reported]

colnames(exp_prof_tumor) = df.stages$short.id
colnames(exp_prof_tumor_reported) = df.stages$short.id[indexes.stages.reported]


dds_tumor = DESeqDataSetFromMatrix(countData = exp_prof_tumor, colData = df.stages, design = ~stage)
dds_tumor_reported = DESeqDataSetFromMatrix(countData = exp_prof_tumor_reported, 
                                            colData = df.stages[indexes.stages.reported,], design = ~stage)
dds_tumor = dds_tumor[match(rownames(dds), rownames(dds_tumor)),]

dds_tumor <- dds_tumor[ rowSums(counts(dds_tumor)) > 1, ]
dds_tumor_reported <- dds_tumor_reported[ rowSums(counts(dds_tumor_reported)) > 1, ]

dds_tumor <- DESeq(dds_tumor)
dds_tumor_reported <- DESeq(dds_tumor_reported)

save(dds_tumor, file = 'environment/dds_tumor.RData')
save(dds_tumor_reported, file = 'environment/dds_tumor_reported.RData')


exp_fpqm_tumor <- exp_fpqm[match(remove.dots(rownames(dds_tumor)), rownames(exp_fpqm)), tumor.indexes]
exp_fpqm_tumor_reported <- exp_fpqm_tumor[match(rownames(dds_tumor_reported), rownames(exp_fpqm_tumor)),indexes.stages.reported]

colnames(exp_fpqm_tumor_reported) <- df.stages$short.id[indexes.stages.reported]
exp_fpqm_tumor_log_reported <- log2(exp_fpqm_tumor_reported+1)

exp_prof_tumor_reported = exp_prof_tumor[,indexes.stages.reported]

only.tumor.reported = list() ##Contains the various normalisations of tumor samples for which stages are known
##Also contains info about genes with various categories
only.tumor.reported[['dfs']] = list()
only.tumor.reported[['dfs']][['vs']] = vst(dds_tumor_reported)
only.tumor.reported[['dfs']][['nt']] = normTransform(dds_tumor_reported)
only.tumor.reported[['dfs']][['fpqm']] = as.matrix(exp_fpqm_tumor_reported)
only.tumor.reported[['dfs']][['fpqm_log']] = as.matrix(log2(exp_fpqm_tumor_reported+1))

only.tumor.reported[['genes']] = list()
only.tumor.reported[['genes']] = create.list.imp.genes(only.tumor.reported[['dfs']], 5000)

only.tumor = list()
only.tumor[['dfs']] = list()
only.tumor[['dfs']][['vs']] = vst(dds_tumor)
only.tumor[['dfs']][['nt']] = normTransform(dds_tumor)
only.tumor[['dfs']][['fpqm']] = as.matrix(exp_fpqm_tumor)
only.tumor[['dfs']][['fpqm_log']] = as.matrix(log2(exp_fpqm_tumor + 1))

only.tumor[['genes']] = list()
only.tumor[['genes']] = create.list.imp.genes(only.tumor[['dfs']], 5000)


pca_tum_rep = list() ##stores the pca list for tumor reported for various dfs for various genes
pca_tum_rep = get.list.plotPCA(only.tumor.reported$dfs, only.tumor.reported$genes, 
                               colData = df.stage.tumor.rep, intgroup = 'stage')
pca_tum_rep[['diff']] = get.list.plotPCA.gene(c(2:5), only.tumor.reported$dfs , diff.genes, 'diff', type = 1, 
                                              intgroup = 'stage', colData = df.stage.tumor.rep)

##Folder where images will be saved
image.direct.main = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/images/tumor/PCA'
##Saves the various images
save.plots.under.main.wd(image.direct.main, names(pca_tum_rep), pca_tum_rep, c(3,3,3,2))

lda_tum_rep = list()
lda_tum_rep = get.list.lda(only.tumor.reported$dfs, only.tumor.reported$genes, 
                                    groups = df.stage.tumor.rep$stage)
lda_tum_rep_3 = get.list.lda(only.tumor.reported$dfs, only.tumor.reported$genes, 
                           groups = df.stage.tumor.rep$stage)

lda_tum_rep_3$diff = get.list.lda.gene(c(2:5), only.tumor.reported$dfs, diff.genes, 'diff', type = 1, df.stage.tumor.rep$stage)  
image.direct.main = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/images/tumor/LDA3'
save.plots.under.main.wd(image.direct.main, names(pca_tum_rep), lda_tum_rep, c(3,3,3,2))

entire.data.lda <- lapply(only.tumor.reported$dfs, function(x)
    {
    if(typeof(x) == 'S4')
      x = assay(x)
    x = t(x)
    lda1 = lda(x, grouping = df.stage.tumor.rep$stage)
    lda1.p = predict(lda1)
    lda1.p$x = data.frame(lda1.p$x)
    lda1.p$x$group = df.stage.tumor.rep$stage
    lda1.list = list(lda1.p, lda1)
  })

##lda1 = lda(x = t(assay(only.tumor.reported$dfs$vs)), grouping = df.stage.tumor.rep$stage)


##lda1.p = predict(lda1)
##lda_tum_rep$`100`$x = data.frame(lda_tum_rep$`100`$x)
##lda_tum_rep$`100`$x$group = df.stage.tumor.rep$stage
##View(lda1.p$x)
##View(lda1$means)
##ggplot(data = lda_tum_rep$`100`$x, aes_string(x='LD1', y='LD2', color ='group')) +
##coord_fixed() + geom_point(size=2.5) 
##Test
# library(FactoMineR)
# library(factoextra)
# df1 = t(only.tumor.reported$dfs$fpqm[only.tumor.reported$genes$fpqm$exp[1:100], ])
# df2 = t(only.tumor.reported$dfs$fpqm[only.tumor.reported$genes$fpqm$exp[1:1000], ])
# rownames(df1) = df.stage.tumor.rep$stage
# rownames(df2) = df.stage.tumor.rep$stage
# pc = PCA(df1)
# pc2 = PCA(df2)
# fviz_pca_ind(pc)
# fviz_pca_ind(pc2)
