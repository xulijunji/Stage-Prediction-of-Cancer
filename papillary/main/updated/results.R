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
###############DeSeq2 feature