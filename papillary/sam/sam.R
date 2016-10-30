install.packages('samr')
library(samr)

y.type = c()
j = 1
for(i in seq_along(sample.info$type))
{
  if(i %% 2 == 1)
    y.type = c(y.type,-j)
  else
  {
    y.type = c(y.type,j)
    j = j + 1
  }
}
remove(j)
dds = dds[-which(rowSums(assay(dds)) < 2),]  
res.sam <- SAMseq(x=assay(dds), y = y.type, resp.type = 'Two class paired' )


for(i in seq(2))
{
  res.sam[[4]][[i]] = data.frame(GeneId = res.sam[[4]][[i]][,1],
                                      GeneName = res.sam[[4]][[i]][,2],
                                      Score = as.numeric(res.sam[[4]][[i]][,3]),
                                      FoldChange = as.numeric(res.sam[[4]][[i]][,4]),
                                      q_value = as.numeric(res.sam[[4]][[i]][,5])
                                    )
}
res.sam <- res.sam.copy

diff.genes.sam <- list()
diff.genes.sam[['up']] = list()
diff.genes.sam[['down']] = list()
diff.genes.sam$up[['0.01_1']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.up[,5] < 1 & 
                                                              res.sam$siggenes.table$genes.up[,4] > 1)]
diff.genes.sam$up[['0.05_1']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.up[,5] < 5 & 
                                                              res.sam$siggenes.table$genes.up[,4] > 1)]
diff.genes.sam$up[['0.05_2']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.up[,5] < 5 & 
                                                              res.sam$siggenes.table$genes.up[,4] > 2)]
diff.genes.sam$up[['0.05_3']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.up[,5] < 5 & 
                                                              res.sam$siggenes.table$genes.up[,4] > 3)]
diff.genes.sam$up[['0.05_4']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.up[,5] < 5 & 
                                                              res.sam$siggenes.table$genes.up[,4] > 4)]
diff.genes.sam$up[['0.05_5']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.up[,5] < 5 & 
                                                              res.sam$siggenes.table$genes.up[,4] > 5)]

diff.genes.sam$down[['0.01_0.9']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.lo[,5] < 1 & 
                                                                  res.sam$siggenes.table$genes.lo[,4] > 0.9)]
diff.genes.sam$down[['0.05_0.9']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.lo[,5] < 5 & 
                                                                  res.sam$siggenes.table$genes.lo[,4] > 0.9)]
diff.genes.sam$down[['0.01_0.99']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.lo[,5] < 1 & 
                                                                  res.sam$siggenes.table$genes.lo[,4] > 0.99)]

diff.genes.sam$down[['0.05_0.95']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.lo[,5] < 5 & 
                                                                   res.sam$siggenes.table$genes.lo[,4] > 0.95)]
diff.genes.sam$down[['0.05_0.99']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.lo[,5] < 5 & 
                                                                   res.sam$siggenes.table$genes.lo[,4] > 0.99)]
diff.genes.sam$down[['0.05_0.97']] =  rownames(assay(dds))[which(res.sam$siggenes.table$genes.lo[,5] < 5 && 
                                                                   res.sam$siggenes.table$genes.lo[,4] > 0.97)]

length(diff.genes.sam$down$`0.05_0.99`)
y.sam <- c()
for(i in stages.levels)
{
  if(i == 'stage i')
    y.sam <- c(y.sam,1)
  else if(i == 'stage ii')
    y.sam <- c(y.sam,2)
  else if(i == 'stage iii')
    y.sam <- c(y.sam,3)
  else
    y.sam <- c(y.sam,4)
}
res.sam.stages <- SAMseq(x=assay(dds), y = stages.levels, resp.type = 'Multiclass' )
