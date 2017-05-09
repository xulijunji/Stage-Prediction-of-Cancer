source('main/updated/initialisation.R')
load('environment/accuracy_feature/updated/cv_model.RData')
load('environment/accuracy_feature/updated/train_model.RData')
load('environment/accuracy_feature/updated/test_pred.RData')
cv.results <- list()
test.results <- list()

##############Shrunken feature#################
########CV########
cv.results[['shrunken']] <- list()

cv.results$shrunken[['shrunken']] <- get.eval.list(stages.levels.comb[train.indexes], 
                                              cv.model$shrunken$shrunken$pred)
knn.pred <- lapply(cv.model$shrunken$knn, function(model)
  {
  model$pred
})
cv.results$shrunken[['knn']] <- get.eval.list(stages.levels.comb[train.indexes], 
                                                   knn.pred)
cv.results$shrunken[c('svm','rf', 'nb')] <- lapply(cv.model$shrunken[c('svm','rf', 'nb')],
                                        function(pred.cv){
                                        get.eval.list(stages.levels.comb[train.indexes], pred.cv)
                                                   })
#########Test######
test.results[['shrunken']] <- list()
test.results$shrunken <- lapply(test.pred$shrunken, function(pred.test){
  get.eval.list(stages.levels.comb[test.indexes], pred.test)
})
#############Shrunken feature##################

###############VarSelRF feature###############
###################CV###########
cv.results[['varSelRF']] <- list()
cv.results$varSelRF[['shrunken']] <- get.eval.list(stages.levels.comb[train.indexes],
                                                   cv.model$varSelRF$shrunken$pred)
knn.pred <- lapply(cv.model$varSelRF$knn, function(model)
{
  model$pred
})
cv.results$varSelRF[['knn']] <- get.eval.list(stages.levels.comb[train.indexes], knn.pred)
cv.results$varSelRF[c('svm','rf', 'nb')] <- lapply(cv.model$varSelRF[c('svm','rf', 'nb')],
                                                   function(pred.cv){
                                                get.eval.list(stages.levels.comb[train.indexes], pred.cv)
                                                   })
##########Test########
test.results[['varSelRF']] <- list()
test.results$varSelRF <- lapply(test.pred$varSelRF, function(pred.test){
  get.eval.list(stages.levels.comb[test.indexes], pred.test)
})
###############VarSelRF feature##############

###############DeSeq2 feature#################
########CV##########
cv.results[['deseq2_2']] <- list()
cv.results$deseq2_2[['shrunken']] <- get.eval.list(stages.levels.comb[train.indexes],
                                                   cv.model$deseq2_2$shrunken$pred)
knn.pred <- lapply(cv.model$deseq2_2$knn, function(model)
{
  model$pred
})
cv.results$deseq2_2[['knn']] <- get.eval.list(stages.levels.comb[train.indexes], knn.pred)
cv.results$deseq2_2[c('svm','rf', 'nb')] <- lapply(cv.model$deseq2_2[c('svm','rf', 'nb')],
                                                   function(pred.cv){
                                           get.eval.list(stages.levels.comb[train.indexes], pred.cv)
                                                   })

cv.results[['deseq2_1.5']] <- list()
cv.results$deseq2_1.5[['shrunken']] <- get.eval.list(stages.levels.comb[train.indexes],
                                                   cv.model$deseq2_1.5$shrunken$pred)
knn.pred <- lapply(cv.model$deseq2_1.5$knn, function(model)
{
  model$pred
})
cv.results$deseq2_1.5[['knn']] <- get.eval.list(stages.levels.comb[train.indexes], knn.pred)
cv.results$deseq2_1.5[c('svm','rf', 'nb')] <- lapply(cv.model$deseq2_1.5[c('svm','rf', 'nb')],
                                                   function(pred.cv){
                                                     get.eval.list(stages.levels.comb[train.indexes], pred.cv)
                                                   })

cv.results[['deseq2_1']] <- list()
cv.results$deseq2_1[['shrunken']] <- get.eval.list(stages.levels.comb[train.indexes],
                                                   cv.model$deseq2_1$shrunken$pred)
knn.pred <- lapply(cv.model$deseq2_1$knn, function(model)
{
  model$pred
})
cv.results$deseq2_1[['knn']] <- get.eval.list(stages.levels.comb[train.indexes], knn.pred)
cv.results$deseq2_1[c('svm','rf', 'nb')] <- lapply(cv.model$deseq2_1[c('svm','rf', 'nb')],
                                                   function(pred.cv){
                                                     get.eval.list(stages.levels.comb[train.indexes], pred.cv)
                                                   })

cv.results$deseq2_2.5[c('svm','rf', 'nb')] <- lapply(cv.model$deseq2_2.5[c('svm','rf', 'nb')],
                                                   function(pred.cv){
                                                     get.eval.list(stages.levels.comb[train.indexes], pred.cv)
                                                   })
cv.results$deseq2_3[c('svm','rf', 'nb')] <- lapply(cv.model$deseq2_3[c('svm','rf', 'nb')],
                                                     function(pred.cv){
                                get.eval.list(stages.levels.comb[train.indexes], pred.cv)
                                                     })

#######Test#########
test.results[['deseq2_2']] <- list()
test.results$deseq2_2 <- lapply(test.pred$deseq2_2, function(pred.test){
  get.eval.list(stages.levels.comb[test.indexes], pred.test)
})
test.results[['deseq2_1.5']] <- list()
test.results$deseq2_1.5 <- lapply(test.pred$deseq2_1.5, function(pred.test){
  get.eval.list(stages.levels.comb[test.indexes], pred.test)
})
test.results[['deseq2_1']] <- list()
test.results$deseq2_1 <- lapply(test.pred$deseq2_1, function(pred.test){
  get.eval.list(stages.levels.comb[test.indexes], pred.test)
})

###############DeSeq2 feature#################
save(cv.results, file = 'environment/accuracy_feature/updated/cv_results.RData')
save(test.results, file = 'environment/accuracy_feature/updated/test_results.RData')

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

  