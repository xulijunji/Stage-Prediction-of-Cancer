###Here I am trying to check the labels which I have so far arrived at actually correlate with
##what is given by the TCGA biolinks

library(TCGAbiolinks)
load('environment/sample_info.RData')
load('environment/sample_info_tumor_rep_normal.RData')
load('environment/dds_tumor_not_match_stage.RData')

q <- GDCquery(project = 'TCGA-KIRP', data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", 
              workflow.type = "HTSeq - Counts")
c <- GDCquery_clinic(project = 'TCGA-KIRP')


get.case.from.sample <- function(samples)
{
  cases <- sapply(strsplit(samples, split = '-', fixed = T), function(case)
  {
    paste(case[1], case[2], case[3], sep = '-')
  })
  return(cases)
}

table(sapply(strsplit(q$results[[1]]$cases, '-'), function(x) x[4]))
##From here we have 289 tumor samples and 32 normal samples

table(c$tumor_stage)
##From here have 291 tumor samples
length(unique(c$submitter_id))

##
library(dplyr)
tumor.rep.cases <- filter(sample.info.all.rep, type == 'T' ) %>% select(sample.names) %>% unlist %>%
 get.case.from.sample

tumor.nor.rep.cases <- filter(sample.info.all.rep, type == 'N') %>% select(sample.names) %>% unlist%>%
  get.case.from.sample

tumor.not.rep.cases.c <- filter(c, tumor_stage == 'not reported') %>% select(submitter_id) %>% unlist
head(tumor.nor.rep.cases)


stagei.cases.c <- filter(c, tumor_stage == 'stage i') %>% select(submitter_id) %>% unlist
stageiii.cases.c <- filter(c, tumor_stage == 'stage iii') %>% select(submitter_id) %>% unlist
stageii.cases.c <- filter(c, tumor_stage == 'stage ii') %>% select(submitter_id) %>% unlist

stagei.cases <- filter(sample.info.all.rep, stage.type == 'stage i') %>% select(sample.names) %>% 
  unlist %>% get.case.from.sample
stageiii.cases <- filter(sample.info.all.rep, stage.type == 'stage iii') %>% select(sample.names) %>% 
  unlist %>% get.case.from.sample
stageii.cases <- filter(sample.info.all.rep, stage.type == 'stage ii') %>% select(sample.names) %>% 
  unlist %>% get.case.from.sample
stageiv.cases <- filter(sample.info.all.rep, stage.type == 'stage iv') %>% select(sample.names) %>% 
  unlist %>% get.case.from.sample

norm.cases

length(intersect(tumor.nor.rep.cases, tumor.not.rep.cases.c)) ###0
length(intersect(tumor.rep.cases, tumor.not.rep.cases.c)) ###0
length(intersect(tumor.rep.cases, tumor.not.rep.cases.c))

q.cases <- get.case.from.sample(q$results[[1]]$cases)
head(q.cases)
length(intersect(q.cases, tumor.rep.cases))  ##259 expected
length(intersect(q.cases, tumor.not.rep.cases.c)) ##29 not 30
length(intersect(q.cases, stagei.cases)) ###172 not 173


##Sorted things work for our HTSeq count data 289 cases not 291 cases
length(intersect(q.cases, stagei.cases)) ##172
length(intersect(stagei.cases.c, stagei.cases)) ##172 out of 173

length(intersect(q.cases, stageiii.cases)) ##51
length(intersect(stageiii.cases.c, stageiii.cases)) ###51
length(intersect(q.cases, stageiii.cases.c)) ###52 not 51
length(intersect(q.cases, stageii.cases))  ##21 not 22 but justified\
length(intersect(q.cases, stageiv.cases)) ##15

setdiff(stageii.cases, stageii.cases.c)

setdiff(stageiii.cases.c, stageiii.cases) %in% tumor.rep.cases   ###False
setdiff(stageiii.cases.c, stageiii.cases) %in% tumor.nor.rep.cases  ###False
##So far for stage iii we have an issue
##Note total count from our set 260 + 31 + 29 + 1 
