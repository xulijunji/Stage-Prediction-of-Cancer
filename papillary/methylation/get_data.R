library(TCGAbiolinks)
load('environment/sample_info_tumor_rep_normal.RData')
query.met <- GDCquery(project = "TCGA-KIRP", data.category = 'DNA Methylation')
GDCdownload(query.met, directory = '../../../temp_kirp/')
data <- GDCprepare(query.met, directory = '../../../temp_kirp')

query.met.450k <- GDCquery(project = c("TCGA-KIRP"),
                                       data.category = "DNA Methylation",
                                       platform = c("Illumina Human Methylation 450"))
  
data.450k <- GDCprepare(query.met.450k, directory = '../../../KIRP/450/')  
table(query.met$results[[1]]$tissue.definition)

head(query.met$access)
table(query.met$results[[1]]$platform)

table(query.met$results[[1]]$platform[tum.inds])
head(query.met$results[[1]]$cases)

sum(query.met$results[[1]]$cases[tum.inds] %in% sample.info.all.rep$sample.names)

library(dplyr)
get.case.from.sample <- function(samples)
{
  cases <- sapply(strsplit(samples, split = '-', fixed = T), function(case)
  {
    paste(case[1], case[2], case[3], sep = '-')
  })
  return(cases)
}

case.tum.rep.rna.ids <- filter(sample.info.all.rep, type == 'T') %>% select(sample.names) %>% unlist
case.tum.rep.rna.ids <- get.case.from.sample(case.tum.rep.rna.ids)
head(case.tum.rep.rna.ids)
length(case.tum.rep.rna.ids) == 260

#case.tum.meth.ids <- filter(query.met$results[[1]], tissue.definition == levels(tissue.definition)[14]) %>%
 # select(cases) %>% unlist
meth.tum.inds <-  which(colData(data)[,4] == 'TP' | colData(data)[,4] == 'TAP')
case.meth.tum.ids <- get.case.from.sample(colnames(data)[meth.tum.inds])
length(case.meth.tum.ids)
head(case.meth.tum.ids)  

sum(case.tum.rep.rna.ids %in% case.meth.tum.ids) == length(case.tum.rep.rna.ids)

meth.rep.tum.inds <- match(case.tum.rep.rna.ids, case.meth.tum.ids)
length(meth.rep.tum.inds)
table(query.met$results[[1]]$platform[meth.rep.tum.inds])
case.meth.rep.tum.ids <- case.meth.tum.ids[meth.rep.tum.inds]
length(intersect(case.meth.rep.tum.ids, case.meth.rep.tum.ids)) ##Two patients same sample but one with 05
setdiff(case.tum.rep.rna.ids, case.meth.ids[meth.rep.tum.inds])
setdiff(case.meth.ids[meth.rep.tum.inds], case.tum.rep.rna.ids)
setdiff(case.meth.rep.tum,case.tum.rep.rna.ids)
which(case.meth.tum.ids == 'TCGA-UZ-A9PS')
colnames(data)[meth.tum.inds][c(140,218)]

reg.int <- filter(sample.info.all.rep, type == 'T')
nrow(reg.int)
sum(droplevels(reg.int$stage.type) == droplevels(stages.levels)) == length(stages.levels)

which(duplicated(meth.rep.tum.inds))
meth.rep.tum.inds[242] <- 218
meth.rep.tum.inds <- meth.tum.inds[meth.rep.tum.inds]

pheno.info.meth.tum.rep <- query.met$results[[1]][meth.rep.tum.inds,]
table(colData(data)[meth.rep.tum.inds,10])
write(pheno.info.meth.tum.rep$file_id, '../../../KIRP/req_files.txt')
table(pheno.info.meth.tum.rep$platform)

meth.27.ind <- which(pheno.info.meth.tum.rep$platform == levels(as.factor(pheno.info.meth.tum.rep$platform))[1])
stages.levels.comb[meth27.array]

save(data, file = 'environment/methylation/whole_data.RData')
save(pheno.info.meth.tum.rep, file = 'environment/methylation/pheno_meth_tum_rep.RData')
save(meth.rep.tum.inds, file = 'environment/methylation/meth_tum_rep_inds.RData')
