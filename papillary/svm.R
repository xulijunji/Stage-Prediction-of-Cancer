####testing file to be forgotten, using it for learning purposes
install.packages('e1071')
library(e1071)
library(caret)
exp_fpqm_tumor_reported_diff_5 <- as.matrix(t.data.frame(exp_fpqm_tumor_reported[diff.genes[[5]], ]))
svm.whole <- svm(x = exp_fpqm_tumor_reported_diff_5, y = stages.levels, kernel = 'radial', gamma = 1)
pred.svm = predict(svm.whole, exp_fpqm_tumor_reported_diff_5)
table(stages.levels, pred.svm)


cv.svm(exp_fpqm_tumor_reported_diff_5, 10)

library(class)
pred = knn(train = exp_fpqm_tumor_reported_diff_5, test = exp_fpqm_tumor_reported_diff_5, cl = stages.levels, k = 1)
table(stages.levels, pred)



remove(exp_fpqm_tumor_reported_diff_5)
remove(pred.svm, pred,svm_tune)
pred.knn <- knn.cv(exp_fpqm_tumor_reported_diff_5, cl = stages.levels, k = 5)
table(stages.levels, pred.knn)

svm_tune <- tune(svm, train.x=exp_fpqm_tumor_reported_diff_5, train.y=stages.levels, 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))