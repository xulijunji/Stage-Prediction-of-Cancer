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
diff = maxs-mins
hist(diff)
hist(sds)

max.counts <- apply(assay(dds_tumor_reported), 2, max)
min.counts <- apply(assay(dds_tumor_reported), 2, min)
diff.counts <- max.counts - min.counts

maxs.vs <- apply(req.dfs$vs, 2, max)
mins.vs <- apply(req.dfs$vs, 2, min)
diff.vs <- maxs.vs - mins.vs
sds.vs <- apply(req.dfs$vs, 2, sd)
hist(diff.vs)
hist(sds.vs)
sum(sds.vs < 0.3)
sum(diff.vs < 2)
length(intersect(which(diff.vs < 2), which(sds.vs < 0.3)))

genes.to.removed.ids.sds <- which(sds < 1e4)
sum(sds < 1e4)
length(intersect(which(diff < 1e5), which(sds < 1e4)))
length(intersect(which(sds < 1e4), stable_genes$fpqm$cho$V1[1:150]))

req.dfs.rem = req.dfs
req.dfs.rem$fpqm = req.dfs$fpqm[,-intersect(which(sds < 1e4), which(diff < 1e5))]
req.dfs.rem$vs = req.dfs.rem$vs[,-intersect(which(diff.vs < 2), which(sds.vs < 0.3))]
wd = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/data/'
write.csv(req.dfs.rem$vs, paste(wd, 'vs_proc.csv', sep = ''))
write.csv(req.dfs.rem$fpqm, paste(wd, 'fpqm_proc.csv', sep = ''))

View(req.dfs$fpqm[,order(mins)[1:9750]])

stable_genes = list()
wd = '~/Dropbox/honours/sem 7/RNA_Seq/papillary/results/tumor/stable_gene_sel/'
stable_genes_files = list.files(wd)
names.stable.genes <- c('1','2','cho','f')
names.dfs <- c('fpqm','fpqm_log','nt','vs')
stable_genes_files = stable_genes_files[2:17]
j = 1
for(i in stable_genes_files)
{
  quo = ceiling(j / 4)
  rem = j %% 4
  if(rem == 0)
  {
    rem = 4
    f = read.csv(paste(wd,i,sep = ''), header = T)
  }
  else
    f = read.csv(paste(wd,i,sep = ''), header = F)
  
  stable_genes[[names.dfs[quo]]][[names.stable.genes[rem]]] = f
  j = j + 1
}
stable_genes$fpqm$`1`$V1
