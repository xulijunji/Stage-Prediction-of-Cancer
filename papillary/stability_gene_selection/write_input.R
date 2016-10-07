load('~/Dropbox/honours/sem 7/RNA_Seq/papillary/environment/only_tumor_reported.RData')
load('~/Dropbox/honours/sem 7/RNA_Seq/papillary/environment/dds_tumor_reported.RData')
req.dfs <- lapply(only.tumor.reported$dfs, function(x)
{
  library(DESeq2)
  if(typeof(x) == 'S4')
    t.data.frame(as.data.frame(assay(x)))
  else
    t.data.frame(as.data.frame(x))
})
names(req.dfs) <- names(only.tumor.reported$dfs)
write.dfs.csvs('data',req.dfs)

View(req.dfs$fpqm)
mins <- apply(req.dfs$fpqm, 2, min)
sort(mins)[1:9760]
mins[order(mins)[1:9750]]
maxs<- apply(req.dfs$fpqm, 2, max)
vars <- apply(req.dfs$fpqm, 2, var)
sds <- apply(req.dfs$fpqm, 2, sd)
View(req.dfs$fpqm[,order(mins)[1:9750]])
