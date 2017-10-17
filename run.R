library(rprojroot)
library(tidyverse)
library(stringr)
library(scales)
library(pamr)
library(caret)
library(GEOquery)
library(illuminaHumanv4.db)
library(Matrix)
library(org.Hs.eg.db)

thres = 0 

#grab data
load("predictor.RData")
gse <- getGEO('GSE81423', GSEMatrix = TRUE)
gse <- gse[[1]]
pdata <- pData(gse)
eset.test <- data.frame(exprs(gse))

#grab illumina annotation
anno <- select(illuminaHumanv4.db, keys=rownames(eset.test), 
	keytypes="PROBEID", columns="ENSEMBL")
anno <- anno[!duplicated(anno$PROBEID),]
input.data <-  cbind(eset.test, ensembl=anno$ENSEMBL)
input.data <- input.data[!is.na(input.data$ensembl),]
 
#parse duplicated ensembl id
id.dup <- input.data[duplicated(input.data$ensembl) 
	| duplicated(input.data$ensembl, 
	fromLast = TRUE),ncol(input.data)]
data.dup <- as.matrix(input.data[duplicated(input.data$ensembl) |
	duplicated(input.data$ensembl, fromLast = TRUE),])
ezid.dup <- unique(id.dup)
data.unique <- input.data[!input.data$ensembl %in% id.dup,]
rownames(data.unique) <- data.unique$ensembl
data.unique$ensembl = NULL

#assess CV for multimapped probesets
data.dup2 <- lapply(ezid.dup, function(x) { 
	expr <- data.dup[id.dup == x,]
	filtered <- expr[,1:ncol(input.data)-1]
	class(filtered) <- "numeric"
	if(is.matrix(filtered)){
		sd.values <- apply(filtered, 1, sd)
		mean.values <- apply(filtered, 1, mean)
		CV.values <- sd.values/mean.values
		CV.values <- gsub("NaN","0",CV.values)
		expr <- expr[which(CV.values == max(CV.values))[[1]],]
	} else {
		expr
	}
})

#merged data for CV estimate
merger <- data.frame(do.call(rbind, data.dup2))
rownames(merger) <- merger$ensembl
merger$ensembl = NULL
eset <- data.frame(as.matrix(rbind(data.unique, merger)))

#generate sparse matrix
df <-  data.matrix(t(eset))
sparsed <- Matrix(data=0, nrow=nrow(df), ncol=ncol(train.set)-1,
	dimnames=list(rownames(df), colnames(train.set[,-1])), sparse = T)
sparsed <- data.frame(as.matrix(sparsed))
sparsed[sparsed == 0] <- NA
for(i in colnames(sparsed)){
	try(sparsed[, i] <- df[, i])
}
sparsed <- data.matrix(sparsed)
test.data <- data.frame(label = rep('C2', nrow(sparsed)), sparsed)

## normalize test data using training data
test.set <- predict(preProcVal.trn, test.data)
sum(apply(test.set, 2, function(x) all(is.na(x))))
test.set <- t(test.set[,-1])

## scale the test data into [-1, 1]. Warnings are for missing values and can be ignored
test.set <- apply(test.set, 1, function(x) rescale(x, to = c(-1, 1)))

## replace the missing data with smallest normalized expression value (-1). 
## May not be ideal as the missing gene could simply be due to the absence 
## of the probe on the microarray platform
test.set <- apply(test.set, 2, function(x) ifelse(is.na(x), as.numeric(-1), as.numeric(x)))

# scale
## run prediction
prediction <- as.character(pamr.predict(dat.train, t(test.set), threshold = thres))
pred <- cbind(sampleID = rownames(test.set), predicted_label = prediction )
## obtain prediction probability for each testing sample
prob <- pamr.predict(dat.train, t(test.set), type = 'posterior', threshold = thres)

## predict indeterminate samples, gap is the diff of prob. btw the binary predicted labels
gap <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
prediction.alt <- vector('list', length = length(gap))
for(i in seq_along(gap)){
    prediction.alt[[i]] <- as.character(pamr.indeterminate(prob, mingap = gap[i]))
}
pred.na <- as.data.frame(do.call(cbind, prediction.alt))
colnames(pred.na) <- sapply(gap, function(x) str_c('gap', x, sep = '_'))
rownames(pred.na) <- rownames(test.set)
write.table(pred.na, file="./results.csv", sep="\t", quote=FALSE)
