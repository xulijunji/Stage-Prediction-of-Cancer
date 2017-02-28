source('main/updated/initialisation.R')

cv.results <- list()

############Shrunken feature###############
cv.results[['shrunken']] <- list()

###Shrunken
cv.results$shrunken[['shrunken']] <- shrunken.cv.results(cv.model$shrunken$shrunken, stages.levels.comb,
                                                          train.indexes)
a<-pamr.predict(train.model$shrunken$shrunken$atleast_3,
                as.matrix(t(vst_tumor_tum[test.indexes, net.features.updated$shrunken$atleast_3])), 
                threshold = train.model$shrunken$shrunken$atleast_2$threshold[7])
table(a, stages.levels.comb[test.indexes])
