load('environment/diff_genes.RData')
load('environment/res.RData')
load('environment/dds.RData')
load('environment/sample_info.RData')
source('main/updated/initialisation.R')
source('../Feature_Extract.R')
rownames(res) <- remove.dots(rownames(res))

get.same.sign.genes <- function(res.df, genes.list, sign) ##get either all positive or negative fold genes
{
  library(dplyr)
  res.df <- data.frame(res.df)
  genes <- lapply(genes.list, function(genes)
  {
    if(sign == '+')
      intersect(rownames(res.df[res.df[, 2] > 0, ]), genes)
    else
      intersect(rownames(res.df[res.df[, 2] < 0, ]), genes)
  })
  return(genes)
}
rownames(dds) <- remove.dots(rownames(dds))
##Normal vs Tumor
diff.genes.up <- get.same.sign.genes(res, diff.genes, '+')
diff.genes.low <- get.same.sign.genes(res, diff.genes, '-')

create.boxplots(diff.genes.up$`5`[1], t(assay(vst(dds))), sample.info$type)
sapply(diff.genes, length) ==
mapply(function(x,y) sum(length(x), length(y)), diff.genes.up, diff.genes.low )

sapply(diff.genes.low, length)
sapply(diff.genes.up, length)

#### Within tumor ######
load('environment/accuracy_feature/updated/vst_tum_rep.RData')
load('environment/stages.level.comb.RData')
load('environment/accuracy_feature/updated/new_data/fea.RData')
tumor_diff <- colMeans(vst_tumor_tum[stages.levels.comb == 'stage i',]) - 
  colMeans(vst_tumor_tum[stages.levels.comb == 'stage iv',])
sum(tumor_diff < 0)

genes.up <- names(which(tumor_diff < 0))
genes.low <- names(which(tumor_diff > 0))

fea.trial.list.up <- lapply(fea.trial.list, function(genes.list)
  {
  sapply(genes.list, function(genes) intersect(genes, genes.up))
})

fea.trial.list.low <- lapply(fea.trial.list, function(genes.list)
{
  sapply(genes.list, function(genes) intersect(genes, genes.low))
})

fea.at3.up <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 3))
fea.at4.up <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 4))
fea.at1.up <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 1))
fea.at2.up <- Reduce(intersect, create.list.venn(fea.trial.list.up, 1, 2))

fea.at3.un.up <- Reduce(union, create.list.venn(fea.trial.list.up, 1, 3))
fea.at4.un.up <- Reduce(union, create.list.venn(fea.trial.list.up, 1, 4))

fea.at3.un.low <- Reduce(union, create.list.venn(fea.trial.list.low, 1, 3))
fea.at4.un.low <- Reduce(union, create.list.venn(fea.trial.list.low, 1, 4))

length(intersect(fea.at3.up, diff.genes.low$`2`))

length(intersect(fea.at4.un.low, diff.genes.up$`2`))
length(intersect(fea.at3.un.low, diff.genes.up$`2`))

load('environment/genes_map.RData')
write.genes <- function(genes.list, genes.map, path)
{
  lapply(seq_along(genes.list), function(i) 
  {
    write.table(genes.map[genes.list[[i]],], 
              file = paste(path, paste(names(genes.list)[i], 'csv', sep = '.'), sep = '/'), 
              sep = '\t', row.names = F)    
  })
}
write.genes(diff.genes.up, genes.map, 'results/tumorvsnormal/up')
write.genes(diff.genes.low, genes.map, 'results/tumorvsnormal/low')
write.genes(list('at3_up' = fea.at3.up, 'at4_up' = fea.at4.up), genes.map, 
            'results/tumor/int')
write.genes(list('at3_up' = fea.at3.un.up, 'at4_up' = fea.at4.un.up, 
                 'at3_low' = fea.at3.un.low, 'at4_low' = fea.at4.un.low), genes.map, 
            'results/tumor/union')
write.genes(list('at1_up' = fea.at1.up, 'at2_up' = fea.at2.up), genes.map, 
            'results/tumor/int')
write.genes(list(intersect(diff.genes$`2`, Reduce(intersect, d))))