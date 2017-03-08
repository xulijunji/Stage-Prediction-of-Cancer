library(hgu133plus2.db)
library(matrixStats)
dexp1 <- exprs(d)
dexp2 <- as.data.frame(d)
#rownames(dexp2) <- rownames(dexp)
dexp2$ID <- rownames(dexp2)
a <- as.data.frame(hgu133plus2ENTREZID)
colnames(hab) <- c("ID","Entrez")


freq <- as.data.frame(table(a$ID))
freq <- freq[(freq$Freq<2),]
a <- a[a$ID %in% freq$Var1,]
rownames(a) <- a$ID
a <- na.omit(a)

#b <- as.data.frame(table(a$ID))
dexp3 <- dexp2[hab$ID,]
dexp3 <- merge(hab, dexp3, by='ID')

dexp3$IQR <- rowIQRs(dexp3[,c(3:31)])
dexp3$SD <- rowSds(as.matrix(dexp3[,c(1:65)]))
dexp3$MEAN <- rowMeans(as.matrix(dexp3[,c(1:65)]))
dexp3$CV <- dexp3$SD/dexp3$MEAN


#dexp31$CV <- apply(t(dexp31[,2:32]),1,co.var)
dexp4 <- dexp3[order(dexp3$Entrez,dexp3$IQR, decreasing = TRUE),]
#dexp5 <- aggregate(. ~ Entrez, data=dexp4[2:32], FUN=mean)
dexp5 <- dexp4[!duplicated(dexp4$Entrez),]
write.table(dexp5, "C-AD_a_b_1.02.17.csv", quote = FALSE, sep = ",", row.names = FALSE)