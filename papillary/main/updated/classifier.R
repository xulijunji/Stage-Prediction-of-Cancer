source('main/updated/initialisation.R')

get.train.model <- function(data, tr.ind, fea.list, stages)
{
  train.model <- list()  
  fea.list <- get.class.fea(net.fea)
  for(i in seq_along(fea.list))
  {
    fea.name <- names(fea.list)[i]
    if(is.null(train.model[[fea.name]]))
          train.model[[fea.name]] <- list()
    
    train.model[[fea.name]][['shrunken']] <- build.shrunken.classifier(data, tr.ind, fea.list[[fea.name]],
                                                                         stages)
    train.model[[fea.name]][['rf']] <- build.rf.classifier(data, tr.ind, fea.list[[fea.name]],
                                                           stages)
    train.model[[fea.name]][['nb']] <- build.nb.classifier(data, tr.ind, fea.list[[fea.name]],
                                                           stages)
    train.model[[fea.name]][['svm']] <- build.svm.classifier(data, tr.ind, fea.list[[fea.name]],
                                                           stages)
  }
  return(train.model)
}

get.cv.model <- function(data, tr.ind, fea.list, stages, tr.model)
{
  cv.model <- list()
  fea.list <- get.class.fea(net.fea)
  for(i in seq_along(fea.list))
  {
    fea.name <- names(fea.list)[i]
    if(is.null(cv.model[[fea.name]]))
      cv.model[[fea.name]] <- list()
    
    cv.model[[fea.name]][['shrunken']] <- cv.shrunken(data, 10, fea.list[[fea.name]],
                                                      tr.model[[fea.name]]$shrunken, tr.ind, 
                                                      stages)
    cv.model[[fea.name]][['rf']] <- cv.rf.list(data, 10, fea.list[[fea.name]], tr.ind, stages)
    cv.model[[fea.name]][['nb']] <- cv.nb.list(data, 10, fea.list[[fea.name]], tr.ind, stages)
    cv.model[[fea.name]][['svm']] <- cv.svm.list(data, 10, fea.list[[fea.name]], tr.ind, stages)
    cv.model[[fea.name]][['knn']] <- cv.knn.list(data, 10, fea.list[[fea.name]], tr.ind, stages)
  }
  return(cv.model)
}

get.test.pred <- function(data, tr.ind, te.ind, fea.list, stages, tr.model, cv.model)
{
  test.pred <- list()
  for(i in seq_along(fea.list))
  {
    fea.name <- names(fea.list)[i]
    if(is.null(test.pred[[fea.name]]))
      test.pred[[fea.name]] <- list()
    
    test.pred[[fea.name]][['shrunken']] <- predict.shrunken(tr.model[[fea.name]][['shrunken']],
                                                            fea.list[[fea.name]], data[te.ind, ],
                                                            cv.model[[fea.name]][['shrunken']]$thr)
    test.pred[[fea.name]][['rf']] <- predict.model(tr.model[[fea.name]][['rf']], fea.list[[fea.name]],
                                                      data[te.ind, ])
    test.pred[[fea.name]][['nb']] <- predict.model(tr.model[[fea.name]][['nb']], fea.list[[fea.name]],
                                                   data[te.ind,])
    test.pred[[fea.name]][['svm']] <- predict.model(tr.model[[fea.name]][['svm']], fea.list[[fea.name]],
                                                    data[te.ind,])
    test.pred[[fea.name]][['knn']] <- predict.knn(cv.model[[fea.name]][['knn']], fea.list[[fea.name]],
                                                 data[tr.ind,], data[te.ind,], stages[tr.ind])
  }
  return(test.pred)
}

fea.trial.list <-  get.class.fea(net.features.trial)
train.trial.model <- get.train.model(vst_tumor_tum, train.trial.ind, fea.list, stages.levels.comb)
cv.trial.model <- get.cv.model(vst_tumor_tum, train.trial.ind, fea.list, stages.levels.comb, tr.model)
test.pred.trial <- get.test.pred(vst_tumor_tum, train.trial.ind, test.trial.ind, fea.list, 
                                 stages.levels.comb, train.trial.model, cv.trial.model)

save(train.trial.model, file = 'environment/accuracy_feature/updated/new_data/train_model.RData')
save(cv.trial.model, file = 'environment/accuracy_feature/updated/new_data/cv_model.RData')
save(test.pred.trial, file = 'environment/accuracy_feature/updated/new_data/test_pred.RData')
