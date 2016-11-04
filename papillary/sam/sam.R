load('environment/dds.RData')
load('environment/dds_tumor_reported.RData')
load('environment/res.RData')
load('environment/diff_genes.RData')
load('environment/res_sam.RData')
load('environment/res_sam_stages.RData')
load('environment/only_tumor_reported.RData')
load('environment/stages_levels.RData')
load('environment/res_sam_pair_wise.RData')
source('../Feature_Extract.R')
source('sam/sam_func_cop.R')
source('sam/samr.morefuns.R')

install.packages('samr')
library(samr)
library(VennDiagram)

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
                                      q_value = as.numeric(res.sam[[4]][[i]][,5]),
                                      diff = log(as.numeric(res.sam[[4]][[i]][,4]), base = 2)
                                    )
}

res.sam$siggenes.table$comb = rbind(res.sam$siggenes.table$genes.up, res.sam$siggenes.table$genes.lo)
diff.genes.sam <- list()
View(assay(dds)[res.sam$siggenes.table$comb$GeneName[res.sam$siggenes.table$comb$log2FC == -Inf],])
get.genes <- function(res, diff, q)
{
  genes = res[,2][abs(res$diff) > diff & res$q_val < q]
  return(as.numeric(levels(genes))[genes])
}

diff.genes.sam[['0.01_2']] =  get.genes(res.sam$siggenes.table$comb, 2, 1)
diff.genes.sam[['0.05_2']] =  get.genes(res.sam$siggenes.table$comb, 2, 5)
diff.genes.sam[['0.01_4']] =  get.genes(res.sam$siggenes.table$comb, 4, 1)
diff.genes.sam[['0.01_8']] =  get.genes(res.sam$siggenes.table$comb, 8, 1)
diff.genes.sam[['0.01_16']] = get.genes(res.sam$siggenes.table$comb, 16, 1) 
diff.genes.sam[['0.01_32']] = get.genes(res.sam$siggenes.table$comb, 32, 1) 
diff.genes.sam[['0.01_64']] = get.genes(res.sam$siggenes.table$comb, 64, 1) 
########Across Stages########################
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

dds_tumor_reported <- dds_tumor_reported[rowSums(assay(dds_tumor_reported)) > 2,]
res.sam.stages <- SAMseq(x=assay(dds_tumor_reported), y = y.sam, resp.type = 'Multiclass' )

res.sam.stages$siggenes.table$genes.up <- data.frame(GeneID = res.sam.stages$siggenes.table$genes.up[,1],
                                                     GeneName = res.sam.stages$siggenes.table$genes.up[,2],
                                                     Score = res.sam.stages$siggenes.table$genes.up[,3],
                                                     cont1 = as.numeric(res.sam.stages$siggenes.table$genes.up[,4]),
                                                     cont2 = as.numeric(res.sam.stages$siggenes.table$genes.up[,5]),
                                                     cont3 = as.numeric(res.sam.stages$siggenes.table$genes.up[,6]),
                                                     cont4 = as.numeric(res.sam.stages$siggenes.table$genes.up[,7]),
                                                     qvalue = as.numeric(res.sam.stages$siggenes.table$genes.up[,8]))
diff.genes.sam.stages <- list()
diff.genes.sam.stages[[1]] = res.sam.stages$siggenes.table$genes.up$GeneName[
                                    abs(res.sam.stages$siggenes.table$genes.up$cont1) > 7 &
                                      res.sam.stages$siggenes.table$genes.up$qvalue < 1]
diff.genes.sam.stages[[2]] = res.sam.stages$siggenes.table$genes.up$GeneName[
                                  abs(res.sam.stages$siggenes.table$genes.up$cont2) > 3 &
                                  res.sam.stages$siggenes.table$genes.up$qvalue < 1]

diff.genes.sam.stages[[3]] = res.sam.stages$siggenes.table$genes.up$GeneName[
  abs(res.sam.stages$siggenes.table$genes.up$cont3) > 6 &
    res.sam.stages$siggenes.table$genes.up$qvalue < 1]

diff.genes.sam.stages[[4]] = res.sam.stages$siggenes.table$genes.up$GeneName[
  abs(res.sam.stages$siggenes.table$genes.up$cont4) > 4.5 &
    res.sam.stages$siggenes.table$genes.up$qvalue < 1]

diff.genes.sam.stages = lapply(diff.genes.sam.stages, function(x)
  {
  as.numeric(levels(x))[x]
})

names(diff.genes.sam.stages) = c(1,2,3,4)

venn.diagram(diff.genes.sam.stages, filename = 'images/tumor/sam/7_3_6_4.5.tiff')

stage.ind <- sapply(levels(stages.levels), function(x)
                    {which(stages.levels == x)})

get.sam.pairwise <- function(data, main, to.comp, y.sam, stage.ind)
{
  pair.wise <- lapply(to.comp, function(x)
  {
  indexes = Reduce(union,stage.ind[c(main,x)])
  y = y.sam[indexes]
  main.index = which(y == main)
  x.index = which(y == x)
  y[main.index] = 1
  y[x.index] = 2
  print(y)
  SAMseq(x=data[,indexes], y = y, resp.type = 'Two class unpaired') 
  })
  names(pair.wise) = to.comp
  return(pair.wise)
}

get.comb.sam <- function(res.sam)
{
  res.sam.pair.wise.comb <- lapply(res.sam, function(x)
  {
    req.names = names(x[[4]])
    res.sam.comb = NULL
    flag.up = 0
    flag.lo = 1
    if('genes.up'  %in% req.names)
      flag.up = 1
    if('genes.lo'  %in% req.names)
      flag.lo = 1
    if(flag.lo & flag.up)
      res.sam.comb = rbind(x[[4]][['genes.up']], x[[4]][['genes.lo']])
    else if(flag.lo)
      res.sam.comb = x[[4]][['genes.lo']]
    else if(flag.up)
      res.sam.comb = x[[4]][['genes.up']]
  })
  names(res.sam.pair.wise.comb) = names(res.sam)
  return(res.sam.pair.wise.comb)
}

str.to.num.comb <- function(res.sam)
{
  res.sam.pair.wise.comb <- lapply(res.sam, function(x)
  {
    x = data.frame(GeneId = x[,1],
                   GeneName = as.numeric(x[,2]),
                   Score = as.numeric(x[,3]),
                   FoldChange = as.numeric(x[,4]),
                   q_val = as.numeric(x[,5]),
                   diff = log2(as.numeric(x[,4])))
  })
  names(res.sam.pair.wise.comb) = names(res.sam)
  return(res.sam.pair.wise.comb)
}

res.sam.pair.wise = list()
res.sam.pair.wise[[1]] = get.sam.pairwise(assay(dds_tumor_reported), 1,c(2,3,4), y.sam, stage.ind)
res.sam.pair.wise[[2]] = get.sam.pairwise(assay(dds_tumor_reported), 2,c(3,4), y.sam, stage.ind)
res.sam.pair.wise[[3]] = get.sam.pairwise(assay(dds_tumor_reported), 3,c(4), y.sam, stage.ind)

res.sam.pair.wise.comb = list()
res.sam.pair.wise.comb[[1]] = get.comb.sam(res.sam.pair.wise[[1]])
res.sam.pair.wise.comb[[2]] = get.comb.sam(res.sam.pair.wise[[2]])
res.sam.pair.wise.comb[[3]] = get.comb.sam(res.sam.pair.wise[[3]])
remove(res.sam.pair.wise)

res.sam.pair.wise.comb[[1]] = str.to.num.comb(res.sam.pair.wise.comb[[1]])
res.sam.pair.wise.comb[[2]] = str.to.num.comb(res.sam.pair.wise.comb[[2]])
res.sam.pair.wise.comb[[3]] = str.to.num.comb(res.sam.pair.wise.comb[[3]])



genes.list.sam = list()
genes.list.sam[[1]] = list()
genes.list.sam[[1]][['stageii']] = get.genes(res.sam.pair.wise.comb[[1]]$`2`, 1, 1)
genes.list.sam[[1]][['stageiii']] = get.genes(res.sam.pair.wise.comb[[1]]$`3`, 32, 1)
genes.list.sam[[1]][['stageiv']] = get.genes(res.sam.pair.wise.comb[[1]]$`4`, 32, 1)

venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_1_1.tiff')
venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_2_2.tiff')
venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_3_3.tiff')
venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_4_4.tiff')
venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_5_5.tiff')
venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_8_8.tiff')
venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_16_16.tiff')
venn.diagram(genes.list.sam[[1]], filename = 'images/tumor/sam/stages/stagei/1_32_32.tiff')

genes.list.sam[[2]] = list()
genes.list.sam[[2]][['stagei']] = get.genes(res.sam.pair.wise.comb[[1]]$`2`, 1, 1)
genes.list.sam[[2]][['stageiii']] = get.genes(res.sam.pair.wise.comb[[2]]$`3`, 4, 5)
genes.list.sam[[2]][['stageiv']] = get.genes(res.sam.pair.wise.comb[[2]]$`4`, 4, 5)
venn.diagram(genes.list.sam[[2]], filename = 'images/tumor/sam/stages/stageii/1_1_1.tiff')
venn.diagram(genes.list.sam[[2]], filename = 'images/tumor/sam/stages/stageii/1_2_2.tiff')
venn.diagram(genes.list.sam[[2]], filename = 'images/tumor/sam/stages/stageii/1_4_4.tiff')

genes.list.sam[[3]] = list()
genes.list.sam[[3]][['stagei']] = get.genes(res.sam.pair.wise.comb[[1]]$`3`, 8, 1)
genes.list.sam[[3]][['stageii']] = get.genes(res.sam.pair.wise.comb[[2]]$`3`, 2, 1)
genes.list.sam[[3]][['stageiv']] = get.genes(res.sam.pair.wise.comb[[3]]$`4`, 1, 1)
venn.diagram(genes.list.sam[[3]], filename = 'images/tumor/sam/stages/stageiii/8_2_1.tiff')

genes.list.sam[[4]] = list()
genes.list.sam[[4]][['stagei']] = get.genes(res.sam.pair.wise.comb[[1]]$`4`, 8, 1)
genes.list.sam[[4]][['stageii']] = get.genes(res.sam.pair.wise.comb[[2]]$`4`, 1, 1)
genes.list.sam[[4]][['stageiii']] = get.genes(res.sam.pair.wise.comb[[3]]$`4`, 1, 1)
venn.diagram(genes.list.sam[[4]], filename = 'images/tumor/sam/stages/stageiv/1_1_1.tiff')
venn.diagram(genes.list.sam[[4]], filename = 'images/tumor/sam/stages/stageiv/4_1_1.tiff')
venn.diagram(genes.list.sam[[4]], filename = 'images/tumor/sam/stages/stageiv/8_1_1.tiff')

genes.list.sam.
