library(COHCAP)
load('environment/methylation/mat_450k_req.RData')
load('environment/methylation/stage_ind_450k.RData')
source('main/updated/initialisation.R')
####27 k#####
dir.coh = file.path(getwd(),'methylation', 'COHCAP')
project.folder <- file.path(getwd(),'methylation', 'COHCAP')
project.name <- 'cohcap'
write(assay(mat.450k.req), file = file.path(project.folder, 'meth.txt'))
beta.file <- file.path(project.folder, 'meth.csv')
beta.table = COHCAP.annotate(beta.file, project.name, project.folder,
                             platform="450k-UCSC")
filtered.sites <- COHCAP.site(rowData(mat.450k.req), assay(mat.450k.req), p)
write(stages.levels.comb[stage.ind.450k], file = file.path(project.folder, 'stage.txt'), sep = '\t')

filtered.sites = COHCAP.site(file.path(project.folder, 'stage.txt'), rowData(mat.450k.req)[,1:5], project.name,
                             project.folder, ref="parental")
