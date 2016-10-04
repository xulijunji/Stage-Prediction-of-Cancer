####testing file to be forgotten, using it for learning purposes
install.packages('e1071')
library(e1071)

exp_fpqm_tumor_reported_diff_5 <- as.matrix(t.data.frame(exp_fpqm_tumor_reported[diff.genes[[5]], ]))
svm.whole <- svm(x = exp_fpqm_tumor_reported_diff_5, y = stages.levels, kernel = 'radial', gamma = 0.1)
pred.svm = predict(svm.whole, as.matrix(exp_fpqm_tumor_reported_diff_5))
table(stages.levels, pred.svm)

cv.svm <- function(data, folds, gamma)
{
  folds.indexes = createFolds(1:length(rownames(data)), folds)
  #print(folds.indexes)
  for(i in seq_along(folds.indexes))
  {
    train.indexes = unlist(folds.indexes[c(-i)])
    #print(train.indexes)
    test.indexes = unlist(folds.indexes[c(i)])
    length(intersect(train.indexes,test.indexes)) == 0
    svm.whole <- svm(x = data[train.indexes,],
                     y = stages.levels[train.indexes], kernel = 'radial', gamma = gamma)
    pred.svm = predict(svm.whole, data[test.indexes,])
    #print(pred.svm)
    print(table(stages.levels[test.indexes], pred.svm))
  }
}

library(class)
pred = knn(train = exp_fpqm_tumor_reported_diff_5, test = exp_fpqm_tumor_reported_diff_5, cl = stages.levels, k = 1)
table(stages.levels, pred)

library(caret)

remove(exp_fpqm_tumor_reported_diff_5)
pred.knn <- knn.cv(exp_fpqm_tumor_reported_diff_5, cl = stages.levels, k = 5)
table(stages.levels, pred.knn)

svm_tune <- tune(svm, train.x=exp_fpqm_tumor_reported_diff_5, train.y=stages.levels, 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))