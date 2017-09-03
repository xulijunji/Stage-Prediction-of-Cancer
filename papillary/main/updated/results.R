source('main/updated/initialisation.R')
load('environment/accuracy_feature/updated/cv_model.RData')
load('environment/accuracy_feature/updated/train_model.RData')
load('environment/accuracy_feature/updated/test_pred.RData')
cv.results <- list()
test.results <- list()

get.cv.res <- function(tr.ind, cv.model, stages, fea.names)
{
  cv.results <- list()
  for(i in fea.names)
  {
    if(is.null(cv.results[[i]]))
      cv.results[[i]] <- list()
    cv.results[[i]][['shrunken']] <- get.eval.list(stages[tr.ind], 
                                                       cv.model[[i]]$shrunken$pred)
    knn.pred <- lapply(cv.model[[i]][['knn']], function(model)
    {
      model$pred
    })
    cv.results[[i]][['knn']] <- get.eval.list(stages[tr.ind], 
                                                  knn.pred)
    cv.results[[i]][c('svm','rf', 'nb')] <- lapply(cv.model[[i]][c('svm','rf', 'nb')],
                                                       function(pred.cv){
                                                          get.eval.list(stages[tr.ind], pred.cv)
                                                       })
  }   
  return(cv.results)
}

get.test.res <- function(test.pred, te.ind, stages, fea.names)
{
  test.results <- list()
  for(i in fea.names)
  {
    test.results[[i]] <- list()
    test.results[[i]] <- lapply(test.pred[[i]], function(pred.test){
      get.eval.list(stages[te.ind], pred.test)
    })  
  }
  return(test.results)
}
cv.trial.results <- get.cv.res(train.trial.ind, cv.trial.model, stages.levels.comb, names(fea.trial.list))
test.trial.results <- get.test.res(test.pred.trial, test.trial.ind, stages.levels.comb, names(fea.trial.list))
save(cv.trial.results, file = 'environment/accuracy_feature/updated/new_data/cv_results.RData')
save(test.trial.results, file = 'environment/accuracy_feature/updated/new_data/test_results.RData')

###Final Performance
perf.test <- list()
perf.test[['mean']] <- list()
perf.test$mean <- lapply(test.results, function(x)
  {
  publish.results(x, type = 'mean')
})
perf.test$mean$varSelRF <- publish.results(test.results$varSelRF, indexes = c(1,2,3))
perf.test$max <- lapply(test.results, function(x)
{
  publish.results(x, type = 'max')
})
perf.test$max$varSelRF <- publish.results(test.results$varSelRF, indexes = c(1,2,3), type = 'max')

write.dfs.csvs('results/tumor/test/mean/', perf.test$mean)
write.dfs.csvs('results/tumor/test/max/', perf.test$max)

perf.cv <- list()
perf.cv$mean <- lapply(cv.results, function(x)
{
  publish.results(x, type = 'mean')
})
perf.cv$mean$varSelRF <- publish.results(cv.results$varSelRF, indexes = c(1,2,3))

perf.cv$max <- lapply(cv.results, function(x)
{
  publish.results(x, type = 'max')
})
write.dfs.csvs('results/tumor/CV/mean/', perf.cv$mean)
write.dfs.csvs('results/tumor/CV/max/', perf.cv$max)

  