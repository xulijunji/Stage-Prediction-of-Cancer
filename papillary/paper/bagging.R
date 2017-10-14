source('main/updated/initialisation.R')
load('environment/accuracy_feature/updated/new_data/cv_model.RData')
load('environment/accuracy_feature/updated/new_data/test_results.RData')


get.bagging.pred <- function(pred.list)
{
  ##Input pred.list - list of the type containing feature and then prediction with each classifier model
  ##further divided into sets atleast 1,2,3,4
  
  ##Output - of the same form just bagged
  pred.bagged.list <- list()
  for(i in seq_along(pred.list))
  {
    feature.name <- names(pred.list)[i]
    classifiers <- names(pred.list[[feature.name]])
    pred.bagged.list[[feature.name]] <- list()
    for(j in seq_along(pred.list[[feature.name]][[classifiers[1]]]))
    {
      group <- names(pred.list[[feature.name]][[classifiers[1]]])[j]
      stages <- lapply(pred.list[[feature.name]], function(fea)
        {
        fea[[j]]
      })
      
      req.stages <- c()
      for(k in seq_along(stages[[1]]))
      {
        stage.temp <- sapply(stages, function(stage) stage[k])
        counts <- table(stage.temp)
        req.stages <- c(req.stages,names(counts)[which.max(counts)])
      }
      pred.bagged.list[[feature.name]][[group]] <- req.stages
    }
    
  }
  return(pred.bagged.list)
}
a <- get.bagging.pred(test.pred.trial)

table(a$`sam1 fold`$atleast_4, stages.levels.comb[test.trial.ind])
table(test.pred.trial$`sam1 fold`$shrunken$atleast_4  , stages.levels.comb[test.trial.ind])

###End results prelim analysis show no improvement in result rather decrease and expected
###Sens increases and spec decreases
