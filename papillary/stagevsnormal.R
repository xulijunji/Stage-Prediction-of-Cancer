library(DESeq2)
library(VennDiagram)
###Using comp.res from across_tumors.R
load('environment/dds.RData')
load('environment/sample_info.RData')
load('environment/df_stages_rep.RData')
load('environment/exp_prof.RData')
load('environment/diff_genes.RData')
load('environment/stage_comp.RData')
source('../function.R')

match.control.indexes.rep <- match(sample.info$pat.ids, df.stage.tumor.rep$patient.id)
match.control.indexes.rep <- unique(match.control.indexes.rep)

sample.info$type.specifc = as.character('N')
sample.info$type.specifc[sample.info$type == 'T'] <- 
    as.character(df.stage.tumor.rep$stage[match.control.indexes.rep])
sample.info$type.specifc <- as.factor(sample.info$type.specifc)
typeof(sample.info$type.specifc)


exp_prof.diff <- exp_prof[,diff.ids]
colnames(exp_prof.diff) <- rownames(sample.info)
dds_match_control_ind <- DESeqDataSetFromMatrix(exp_prof.diff, colData = sample.info,
                                                design = ~type.specifc)

diff.ids <- match(replace.symbol(sample.info$sample.names, '.', '-'), colnames(exp_prof))

dds_match_control_ind <- dds_match_control_ind[rowSums(assay(dds_match_control_ind)) > 2, ]
dds_match_control_ind <- DESeq(dds_match_control_ind)

results.match.control.stages <- list()
results.match.control.stages[['N']] <- comp.res(dds_match_control_ind,'type.specifc', 'N',
                                                c('stage i','stage ii', 'stage iii', 'stage iv'))
results.match.control.stages[['stage i']] <- comp.res(dds_match_control_ind,'type.specifc', 'stage i',
                                                c('N','stage ii', 'stage iii', 'stage iv'))
results.match.control.stages[['stage ii']] <- comp.res(dds_match_control_ind,'type.specifc', 'stage ii',
                                                       c( 'stage iii', 'stage iv'))
results.match.control.stages[['stage iii']] <- comp.res(dds_match_control_ind,'type.specifc', 'stage iii',
                                                       c( 'stage iv'))

save(results.match.control.stages, file = 'environment/results_match_stages.RData')

genes.list.normal <- list()
genes.list.normal[['normal']] <- list()
genes.list.normal$normal[['stage i']] <- na.omit(get.genes(results.match.control.stages$N$`stage i`,
                                                   7, 0.05, 1))
genes.list.normal$normal[['stage ii']] <- na.omit(get.genes(results.match.control.stages$N$`stage ii`,
                                                           7, 0.05, 1))
genes.list.normal$normal[['stage iii']] <- na.omit(get.genes(results.match.control.stages$N$`stage iii`,
                                                            7, 0.05, 1))
genes.list.normal$normal[['stage iv']] <- na.omit(get.genes(results.match.control.stages$N$`stage iv`,
                                                            7, 0.05, 1))
venn.diagram(genes.list.normal$normal, filename = 'images/tumor/DeSeq2/normalvsstage/normal/7_7_7_7.tiff')

genes.list.normal[['stage i']] <- list()
genes.list.normal$`stage i`[['normal']] <- na.omit(get.genes(results.match.control.stages$N$`stage i`,
                                                             5, 0.05, 1))
genes.list.normal$`stage i`[['stage ii']] <- na.omit(get.genes(results.match.control.stages$`stage i`$`stage ii`,
                                                           4, 0.05, 1))
genes.list.normal$`stage i`[['stage iii']] <- na.omit(get.genes(results.match.control.stages$`stage i`$`stage iii`,
                                                               4, 0.05, 1))
genes.list.normal$`stage i`[['stage iv']] <- na.omit(get.genes(results.match.control.stages$`stage i`$`stage iv`,
                                                               4, 0.05, 1))
venn.diagram(genes.list.normal$`stage i`, filename = 'images/tumor/DeSeq2/normalvsstage/stagei/5_4_4_4.tiff')

genes.list.normal[['stage ii']] <- list()
genes.list.normal$`stage ii`[['normal']] <- na.omit(get.genes(results.match.control.stages$N$`stage ii`,
                                                             5, 0.05, 1))
genes.list.normal$`stage ii`[['stage i']] <- na.omit(get.genes(results.match.control.stages$`stage i`$`stage ii`,
                                                               4, 0.05, 1))
genes.list.normal$`stage ii`[['stage iii']] <- na.omit(get.genes(results.match.control.stages$`stage ii`$`stage iii`,
                                                                4, 0.05, 1))
genes.list.normal$`stage ii`[['stage iv']] <- na.omit(get.genes(results.match.control.stages$`stage ii`$`stage iv`,
                                                               4, 0.05, 1))
venn.diagram(genes.list.normal$`stage ii`, filename = 'images/tumor/DeSeq2/normalvsstage/stageii/5_4_4_4.tiff')

genes.list.normal[['stage iii']] <- list()
genes.list.normal$`stage iii`[['normal']] <- na.omit(get.genes(results.match.control.stages$N$`stage iii`,
                                                              5, 0.05, 1))
genes.list.normal$`stage iii`[['stage i']] <- na.omit(get.genes(results.match.control.stages$`stage i`$`stage iii`,
                                                               5, 0.05, 1))
genes.list.normal$`stage iii`[['stage ii']] <- na.omit(get.genes(results.match.control.stages$`stage ii`$`stage iii`,
                                                                 4, 0.05, 1))
genes.list.normal$`stage iii`[['stage iv']] <- na.omit(get.genes(results.match.control.stages$`stage iii`$`stage iv`,
                                                                4, 0.05, 1))
venn.diagram(genes.list.normal$`stage iii`, filename = 'images/tumor/DeSeq2/normalvsstage/stageiii/5_5_4_4.tiff')

genes.list.normal[['stage iv']] <- list()
genes.list.normal$`stage iv`[['normal']] <- na.omit(get.genes(results.match.control.stages$N$`stage iv`,
                                                               5, 0.05, 1))
genes.list.normal$`stage iv`[['stage i']] <- na.omit(get.genes(results.match.control.stages$`stage i`$`stage iv`,
                                                                4, 0.05, 1))
genes.list.normal$`stage iv`[['stage ii']] <- na.omit(get.genes(results.match.control.stages$`stage ii`$`stage iv`,
                                                                 4, 0.05, 1))
genes.list.normal$`stage iv`[['stage iii']] <- na.omit(get.genes(results.match.control.stages$`stage iii`$`stage iv`,
                                                                 4, 0.05, 1))
venn.diagram(genes.list.normal$`stage iv`, filename = 'images/tumor/DeSeq2/normalvsstage/stageiv/5_4_4_4.tiff')

####On the entire data set excluding not reported
##using info from sample_info.R
indexes.all.rep <- sort(union(tumor.indexes.reported, intersect(normal.indexes, diff.ids)))
sample.info.all.rep <- data.frame(sample.names = sample.ids[indexes.all.rep],
                                  type = type[indexes.all.rep], pat.ids = pats[indexes.all.rep],
                                  stage.type = 'N',
                                  stringsAsFactors = F)
sample.info.all.rep$type = as.factor(sample.info.all.rep$type)

tumor.sample.indexes <- which(sample.info.all.rep$type == 'T')
sample.info.all.rep$stage.type[tumor.sample.indexes] = 
            as.character(df.stages$stage)[match(sample.info.all.rep$sample.names[tumor.sample.indexes], 
                                                       df.stages$sample.id)]
rownames(sample.info.all.rep) = sample.info.all.rep$sample.names
remove(tumor.sample.indexes)
dds_tumor_reported_normal <- DESeqDataSetFromMatrix(countData = exp_prof[,indexes.all.rep],
                                                    colData = sample.info.all.rep, design = ~type)
dds_tumor_reported_normal <- DESeq(dds_tumor_reported_normal)

dds_tumor_reported_normal_stage <- dds_tumor_reported_normal_stage[rowSums(assay(dds_tumor_reported_normal_stage)) > 2, ]
ddds_tumor_reported_normal_stage <- DESeq(dds_tumor_reported_normal_stage)

results.tum.reported.normal <- results(dds_tumor_reported_normal)
diff.genes.tum.rep <- list()
diff.genes.tum.rep[[2]] <- na.omit(get.genes(results.tum.reported.normal, 2, 0.05, 1))
diff.genes.tum.rep[[1]] <- na.omit(get.genes(results.tum.reported.normal, 1, 0.05, 1))
diff.genes.tum.rep[[3]] <- na.omit(get.genes(results.tum.reported.normal, 3, 0.05, 1))
diff.genes.tum.rep[[4]] <- na.omit(get.genes(results.tum.reported.normal, 4, 0.05, 1))
diff.genes.tum.rep[[5]] <- na.omit(get.genes(results.tum.reported.normal, 5, 0.05, 1))


dds_tumor_reported_normal_stage <- DESeqDataSetFromMatrix(countData = exp_prof[,indexes.all.rep],
                                                    colData = sample.info.all.rep, design = ~stage.type)

dds_tumor_reported_normal_stage <- DESeq(dds_tumor_reported_normal_stage)
results.stage.reported.normal <- list()
results.stage.reported.normal[['normal']] <- comp.res(dds_tumor_reported_normal_stage,'stage.type',
                                            'N', c('stage i', 'stage ii', 'stage iii', 'stage iv'))
results.stage.reported.normal[['stage i']] <- comp.res(dds_tumor_reported_normal_stage,'stage.type',
                                                      'stage i', c('stage ii', 'stage iii', 'stage iv'))
results.stage.reported.normal[['stage ii']] <- comp.res(dds_tumor_reported_normal_stage,'stage.type',
                                                       'stage ii', c('stage iii', 'stage iv'))
results.stage.reported.normal[['stage iii']] <- comp.res(dds_tumor_reported_normal_stage,'stage.type',
                                                        'stage iii', c('stage iv'))

genes.list.normal.reported <-list()
genes.list.normal.reported[['N']] <- list()
genes.list.normal.reported$N[['stage i']] <- na.omit(get.genes(results.stage.reported.normal$normal$`stage i`,
                                               5, 0.05, 0.05))
genes.list.normal.reported$N[['stage ii']] <- na.omit(get.genes(results.stage.reported.normal$normal$`stage ii`,
                                                       5, 0.05, 0.05))
genes.list.normal.reported$N[['stage iii']] <- na.omit(get.genes(results.stage.reported.normal$normal$`stage iii`,
                                                       5, 0.05, 0.05))
genes.list.normal.reported$N[['stage iv']] <- na.omit(get.genes(results.stage.reported.normal$normal$`stage iv`,
                                                       5, 0.05, 0.05))
venn.diagram(genes.list.normal.reported$N, 
             filename = 'images/tumor/DeSeq2/reportedvsstage/normal/5_5_5_5.tiff')

genes.list.normal.reported[['stagei']] <- list()
genes.list.normal.reported$stagei[['N']] <- na.omit(get.genes(results.stage.reported.normal$normal$`stage i`,
                                                               2, 0.05, 0.05))
genes.list.normal.reported$stagei[['stage ii']] <- na.omit(get.genes(results.stage.reported.normal$`stage i`$`stage ii`,
                                                                2, 0.05, 0.05))
genes.list.normal.reported$stagei[['stage iii']] <- na.omit(get.genes(results.stage.reported.normal$`stage i`$`stage iii`,
                                                                 2, 0.05, 0.05))
genes.list.normal.reported$stagei[['stage iv']] <- na.omit(get.genes(results.stage.reported.normal$`stage i`$`stage iv`,
                                                                2, 0.05, 0.05))
venn.diagram(genes.list.normal.reported$stagei, 
             filename = 'images/tumor/DeSeq2/reportedvsstage/stagei/2_2_2_2.tiff')

genes.list.normal.reported[['stageii']] <- list()
genes.list.normal.reported$stagei[['N']] <- na.omit(get.genes(results.stage.reported.normal$normal$`stage i`,
                                                              5, 0.05, 0.05))
genes.list.normal.reported$stagei[['stage ii']] <- na.omit(get.genes(results.stage.reported.normal$`stage i`$`stage ii`,
                                                                     3, 0.05, 0.05))
genes.list.normal.reported$stagei[['stage iii']] <- na.omit(get.genes(results.stage.reported.normal$`stage i`$`stage iii`,
                                                                      3.5, 0.05, 0.05))
genes.list.normal.reported$stagei[['stage iv']] <- na.omit(get.genes(results.stage.reported.normal$`stage i`$`stage iv`,
                                                                     3.5, 0.05, 0.05))
venn.diagram(genes.list.normal.reported$stagei, 
             filename = 'images/tumor/DeSeq2/reportedvsstage/stagei/5_3_3.5_3.5.tiff')

###On the tumor not in matched control
indexes.not.match <- sort(union(setdiff(diff.ids,tumor.indexes.reported), 
                                setdiff(tumor.indexes.reported, diff.ids)))
sample.info.not.match <- data.frame(sample.names = sample.ids[indexes.not.match],
                                    type = type[indexes.not.match], pat.ids = pats[indexes.not.match],
                                    stage.type = 'N',
                                    stringsAsFactors = F)
sample.info.not.match$type = as.factor(sample.info.not.match$type)
tumor.indexes.not.match = which(sample.info.not.match$type == 'T')
sample.info.not.match$stage.type[tumor.indexes.not.match] = 
  as.character(df.stages$stage)[match(sample.info.not.match$sample.names[tumor.indexes.not.match], 
                                      df.stages$sample.id)]
sample.info.not.match$stage.type = as.factor(sample.info.not.match$stage.type)
remove(tumor.indexes.not.match)
rownames(sample.info.not.match) = sample.info.not.match$sample.names
dds_tumor_not_match <- DESeqDataSetFromMatrix(countData = exp_prof[,indexes.not.match],
                                              colData = sample.info.not.match, design = ~type)
dds_tumor_not_match_stage <- DESeqDataSetFromMatrix(countData = exp_prof[,indexes.not.match],
                                                    colData = sample.info.not.match, design = ~stage.type)
dds_tumor_not_match <- DESeq(dds_tumor_not_match)
dds_tumor_not_match_stage <- DESeq(dds_tumor_not_match_stage)
dds_tumor_not_match_stage <- dds_tumor_not_match_stage[rowSums(assay(dds_tumor_not_match_stage)) > 2,]
dds_tumor_not_match_stage <- DESeq(dds_tumor_not_match_stage)

results.not.match.stage <- list()
results.not.match.stage[['N']] <- list()
results.not.match.stage[['N']] <- comp.res(dds_tumor_not_match_stage, 'stage.type', 'N',
                                           c('stage i', 'stage ii', 'stage iii', 'stage iv'))
results.not.match.stage[['stage i']] <- comp.res(dds_tumor_not_match_stage, 'stage.type', 'stage i',
                                                 c('stage ii', 'stage iii', 'stage iv'))
results.not.match.stage[['stage ii']] <- comp.res(dds_tumor_not_match_stage, 'stage.type', 'stage ii',
                                                 c('stage iii', 'stage iv'))
results.not.match.stage[['stage iii']] <- comp.res(dds_tumor_not_match_stage, 'stage.type', 'stage iii',
                                                 c('stage iv'))

results.tum.not.match <- results(dds_tumor_not_match)
diff.genes.tum.not.match <- list()
diff.genes.tum.not.match[[2]] <- na.omit(get.genes(results.tum.not.match, 2, 0.05, 1))
diff.genes.tum.not.match[[3]] <- na.omit(get.genes(results.tum.not.match, 3, 0.05, 1))
diff.genes.tum.not.match[[4]] <- na.omit(get.genes(results.tum.not.match, 4, 0.05, 1))
diff.genes.tum.not.match[[5]] <- na.omit(get.genes(results.tum.not.match, 5, 0.05, 1))


genes.list.not.match <- list()
genes.list.not.match[['N']] <- list()
genes.list.not.match$N$'stage i' <- na.omit(get.genes(results.not.match.stage$N$`stage i`, 5, 0.05, 0.05))
genes.list.not.match$N$'stage ii' <- na.omit(get.genes(results.not.match.stage$N$`stage ii`, 5, 0.05, 0.05))
genes.list.not.match$N$'stage iii' <- na.omit(get.genes(results.not.match.stage$N$`stage iii`, 5, 0.05, 0.05))
genes.list.not.match$N$'stage iv' <- na.omit(get.genes(results.not.match.stage$N$`stage iv`, 5, 0.05, 0.05))
venn.diagram(genes.list.not.match$N, filename = 'images/tumor/DeSeq2/not_matchvsstage/normal/5_5_5_5.tiff')

genes.list.not.match[['stage i']] <- list()
genes.list.not.match$`stage i`$'N' <- na.omit(get.genes(results.not.match.stage$N$`stage i`, 5, 0.05, 0.05))
genes.list.not.match$'stage i'$'stage ii' <- na.omit(get.genes(results.not.match.stage$`stage i`$`stage ii`, 3.5, 0.05, 0.05))
genes.list.not.match$'stage i'$'stage iii' <- na.omit(get.genes(results.not.match.stage$`stage i`$`stage iii`, 3, 0.05, 0.05))
genes.list.not.match$`stage i`$'stage iv' <- na.omit(get.genes(results.not.match.stage$`stage i`$`stage iv`, 3.5, 0.05, 0.05))
venn.diagram(genes.list.not.match$`stage i`, filename = 'images/tumor/DeSeq2/not_matchvsstage/stagei/5_2.5_3_3.5.tiff')

genes.list.not.match[['stage ii']] <- list()
genes.list.not.match$`stage ii`$'N' <- na.omit(get.genes(results.not.match.stage$N$`stage ii`, 5, 0.05, 0.05))
genes.list.not.match$'stage ii'$'stage i' <- na.omit(get.genes(results.not.match.stage$`stage i`$`stage ii`, 3, 0.05, 0.05))
genes.list.not.match$'stage ii'$'stage iii' <- na.omit(get.genes(results.not.match.stage$`stage ii`$`stage iii`, 3, 0.05, 0.05))
genes.list.not.match$`stage ii`$'stage iv' <- na.omit(get.genes(results.not.match.stage$`stage ii`$`stage iv`, 3, 0.05, 0.05))
venn.diagram(genes.list.not.match$`stage ii`, filename = 'images/tumor/DeSeq2/not_matchvsstage/stageii/5_3_3_3.tiff')

genes.list.not.match[['stage iii']] <- list()
genes.list.not.match$`stage iii`$'N' <- na.omit(get.genes(results.not.match.stage$N$`stage iii`, 5, 0.05, 0.05))
genes.list.not.match$'stage iii'$'stage i' <- na.omit(get.genes(results.not.match.stage$`stage i`$`stage iii`, 3, 0.05, 0.05))
genes.list.not.match$'stage iii'$'stage ii' <- na.omit(get.genes(results.not.match.stage$`stage ii`$`stage iii`, 3, 0.05, 0.05))
genes.list.not.match$`stage iii`$'stage iv' <- na.omit(get.genes(results.not.match.stage$`stage iii`$`stage iv`, 3, 0.05, 0.05))
venn.diagram(genes.list.not.match$`stage iii`, filename = 'images/tumor/DeSeq2/not_matchvsstage/stageiii/5_3_3_3.tiff')

genes.list.not.match[['stage iv']] <- list()
genes.list.not.match$`stage iv`$'N' <- na.omit(get.genes(results.not.match.stage$N$`stage iv`, 5, 0.05, 0.05))
genes.list.not.match$'stage iv'$'stage i' <- na.omit(get.genes(results.not.match.stage$`stage i`$`stage iv`, 3.5, 0.05, 0.05))
genes.list.not.match$'stage iv'$'stage ii' <- na.omit(get.genes(results.not.match.stage$`stage ii`$`stage iv`, 3, 0.05, 0.05))
genes.list.not.match$`stage iv`$'stage iii' <- na.omit(get.genes(results.not.match.stage$`stage iii`$`stage iv`, 3, 0.05, 0.05))
venn.diagram(genes.list.not.match$`stage iv`, filename = 'images/tumor/DeSeq2/not_matchvsstage/stageiv/5_3.5_3_3.tiff')


venn.diagram(list(diff.genes[[2]], diff.genes.tum.rep[[2]], diff.genes.tum.not.match[[2]]),
             filename = 'images/tumor/DeSeq2/comparing/2_2_2.tiff', 
             category.names = c('matched', 'all_reported', 'not_matched'))
venn.diagram(list(diff.genes[[3]], diff.genes.tum.rep[[3]], diff.genes.tum.not.match[[3]]),
             filename = 'images/tumor/DeSeq2/comparing/3_3_3.tiff', 
             category.names = c('matched', 'all_reported', 'not_matched'))
venn.diagram(list(diff.genes[[4]], diff.genes.tum.rep[[4]], diff.genes.tum.not.match[[4]]),
             filename = 'images/tumor/DeSeq2/comparing/4_4_4.tiff', 
             category.names = c('matched', 'all_reported', 'not_matched'))
venn.diagram(list(diff.genes[[5]], diff.genes.tum.rep[[5]], diff.genes.tum.not.match[[5]]),
             filename = 'images/tumor/DeSeq2/comparing/5_5_5.tiff', 
             category.names = c('matched', 'all_reported', 'not_matched'))
