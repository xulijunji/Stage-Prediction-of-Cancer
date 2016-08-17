df.stage.tumor.rep <- df.stages[indexes.stages.reported, ]
df.stage.tumor.rep$stage = droplevels(df.stage.tumor.rep$stage)
stage.dist = table(df.stage.tumor.rep$stage)

stage.ind = sapply(levels(df.stage.tumor.rep$stage), function(x)
  {
  which(x == df.stage.tumor.rep$stage)
})

max.stage.dist <- max(stage.dist)
total.samples.oversamp <- max.stage.dist*4
indexes.oversamp = (length(df.stage.tumor.rep$stage)+1):total.samples.oversamp

new.indexes.oversamp <- list()

for(i in setdiff(stage.dist, stage.dist[which.max(stage.dist)]))
{
  req.ind <- max.stage.dist - i
  new.indexes.oversamp[[names(stage.dist)[stage.dist == i]]] = sample(indexes.oversamp, req.ind)
  indexes.oversamp <- setdiff(indexes.oversamp, new.indexes.oversamp[[names(stage.dist)[stage.dist == i]]])
}

old.samples.rep <- mapply(function(x,y){
  sample(y,length(x), replace = T)
}, new.indexes.oversamp, stage.ind[-1])

exp_fpqm_tumor_reported_over = t(exp_fpqm_tumor_reported)
exp_fpqm_tumor_reported_over = data.frame(exp_fpqm_tumor_reported_over)

indexes.oversamp = (length(df.stage.tumor.rep$stage)+1):total.samples.oversamp
exp_fpqm_tumor_reported_over[,indexes.oversamp] = rep(1.1, 23731)

mapply(function(x,y)
  {
  #print(x)
  #print(y)
  exp_fpqm_tumor_reported_over[,x] = exp_fpqm_tumor_reported_over[,y]
}, new.indexes.oversamp, old.samples.rep)
for(i in seq_along(new.indexes.oversamp))
{
  for(j in seq_along(old.samples.rep))
  {
    if(i == j)
      exp_fpqm_tumor_reported_over[,new.indexes.oversamp[[i]]] = exp_fpqm_tumor_reported_over[,old.samples.rep[[j]]]
  }
}
exp_fpqm_tumor_reported_over = as.matrix(exp_fpqm_tumor_reported_over)

new.indexes.oversamp$`stage ii`[1:2]
old.samples.rep$`stage ii`[1:2]

stages.oversamp = as.character(df.stage.tumor.rep$stage)
stages.oversamp = c(stages.oversamp, rep('stage i', 428))
stages.oversamp = as.factor(stages.oversamp)
for(i in seq_along(new.indexes.oversamp))             
{
  stages.oversamp[new.indexes.oversamp[[i]]] = names(new.indexes.oversamp)[i]
}
save(exp_fpqm_tumor_reported_over, file = 'environment/exp_fpqm_tumor_reported_over.RData')
save(stages.oversamp, file = 'environment/stages.oversamp')
