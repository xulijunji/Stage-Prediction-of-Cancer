source('main/updated/initialisation.R')
library(FactoMineR)
library(factoextra)



g.set <- net.features.updated$shrunken$atleast_3
p <- do.pca(vst_tumor_tum[train.indexes, g.set], vst_tumor_tum[test.indexes, g.set])
# 
# rf <- randomForest(x = c[[1]], y= stages.levels.comb[train.indexes])
# table(stages.levels.comb[test.indexes],
#       predict(rf, predict(p, vst_tumor_tum[test.indexes, g.set] )[,1:10]))
# rf$confusion
# 
sv <- cv.svm.list(vst_tumor_tum, 10, net.features.updated$shrunken[c(3,4,5,6)], train.indexes, 
                   stages.levels.comb, T)
table(stages.levels.comb[train.indexes], sv[[3]])
sv.train <- svm(x = p[[1]], y = stages.levels.comb[train.indexes])
sv.test <- predict(sv.train, p[[2]])
table(stages.levels.comb[test.indexes], sv.test)

library(kernlab)
kp <- kpca(x = vst_tumor_tum[train.indexes, g.set],
            kernel = 'rbfdot', kpar = list( 0.1))
kp.test <- predict(kp, vst_tumor_tum[test.indexes, g.set])
rf <- randomForest(kp@rotated[,1:10], y = stages.levels.comb[train.indexes])
rf$confusion
pr <- predict(rf, kp.test[,1:10])
table(stages.levels.comb[test.indexes], pr)
View(kp.test)

sv <- ksvm(x = vst_tumor_tum[train.indexes, g.set], y = stages.levels.comb[train.indexes],
           kernel = 'rbfdot', kpar = list(0.01), C = 10)
pr <- predict(sv, vst_tumor_tum[test.indexes, g.set])
table(stages.levels.comb[test.indexes],pr)
