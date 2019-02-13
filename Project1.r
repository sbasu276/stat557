##################################
###### STAT 557 (Project 1) ######
##################################

rm(list=ls()) ## To clear your environment

## Read the data
xTrain=read.csv("ecoli_xTrain.csv", header=FALSE)
yTrain=read.csv("ecoli_yTrain.csv", header=FALSE)
xTest=read.csv("ecoli_xTest.csv", header=FALSE)
yTest=read.csv("ecoli_yTest.csv", header=FALSE)

new_xTrain=read.csv("ecoli_new.xTrain.csv", header=FALSE)
new_yTrain=read.csv("ecoli_new.yTrain.csv", header=FALSE)
new_xTest=read.csv("ecoli_new.xTest.csv", header=FALSE)
new_yTest=read.csv("ecoli_new.yTest.csv", header=FALSE)

#### Part 1 ####
logProd <- function(x){
    lp = sum(x)
    return(lp)
}

logSum <- function(x){
    ls = log(sum(exp(x-max(x))))+max(x)
    return(ls)
}

#### Part 2 ####
prior <- function(yTrain){
    p = as.numeric(table(unlist(yTrain))/nrow(yTrain))
    return(p)
}

likelihood <- function(xTrain, yTrain){
    f = ncol(xTrain) # Number of features
    c = nrow(unique(yTrain)) # Number of classes
    classes = levels(factor(yTrain[,1]))
    df = cbind(xTrain, class=yTrain)
    M = matrix(0, f, c)
    V = matrix(0, f, c)
    colnames(M) = c(classes)
    colnames(V) = c(classes)
    for (j in classes) {
        df1 = df[df[,f+1]==j, 1:f] # Get the observations for class-j
        M[,j] = as.numeric(apply(df1, 2, mean)) # Get mean for each feature
        V[,j] = as.numeric(apply(df1, 2, var)) # Get var for each feature
    }
    return(list(M=M, V=V))
}

naiveBayesClassify <- function(xTest, M, V, p){
    feature = nrow(M)
    cls = ncol(M)
    classes = colnames(M)
    rows = nrow(xTest)
    t = cbind(c(rep(0, rows)))
    for (i in 1:rows) {
        pst = c(rep(0, cls))
        for (c in classes) {
            llh = c(rep(0,feature))
            for (f in 1:feature) {
                llh[f] = dnorm(xTest[i,f], M[f,c], sqrt(V[f,c]), log)
            }
            llh = llh[is.finite(llh)]
            pst[which(classes==c)] = logProd(c(llh, log(p[which(classes==c)])))
        }
        t[i] = as.integer(classes[which.max(pst)])
    }
    return(t)
}

#### Part 3 ####
sigmoidProb <- function(y, x, w){
    d = 1+exp(sum(x*w))
    p = 1/d
    if (y==0) {
        p = 1-p
    }
    return(p)
}

logisticRegressionWeights <- function(xTrain, yTrain, w0, nIter){
    w = w0
    step = 0.1
    f = ncol(xTrain)
    for (itr in 1:nIter) {
        p_hat = apply(xTrain, 1, FUN=function(x) sigmoidProb(1, x, w))
        for (i in 1:f) {
            w[i] = w[i] - step*sum(xTrain[,i]*(yTrain-p_hat))
        }
    }
    return(w)
}

logisticRegressionClassify <- function(xTest, w){
    p = apply(xTest, 1, FUN=function(x) sigmoidProb(1, x, w))
    t = as.integer(p>=0.5)
    t = cbind(t)
    return(t)
}

### Additional Helper functions for Evaluation ###
getPrecision <- function(confMat, class_id) {
    class_id = toString(class_id)
    precision = confMat[class_id, class_id]/sum(confMat[class_id,])
    return(round(precision,digits=3))
}

getRecall <- function(confMat, class_id) {
    class_id = toString(class_id)
    recall = confMat[class_id, class_id]/sum(confMat[,class_id])
    return(round(recall,digits=3))
}

getConfusionMatrix <- function(pred, ref) {
    classes = levels(ref)
    num_cls = length(classes)
    confMat = matrix(0, num_cls, num_cls)
    rownames(confMat) = c(classes)
    colnames(confMat) = c(classes)
    for (i in classes) {
        d = ref[pred==i]
        for (j in classes) {
            confMat[i,j] = sum(d==j)
        }
    }
    return(confMat)
}

evaluate <- function(pred, ref) {
    #pred and ref are two factors for prediction and reference labels
    acc = sum(pred==ref)/length(ref)
    confMat = getConfusionMatrix(pred, ref)
    return(list(accuracy=round(acc, digits=3), confMat=confMat))
}

### Evaluation Part 2 ###
#sink("evaluation.txt", append=FALSE, split=TRUE)
p = prior(yTrain)
L = likelihood(xTrain, yTrain)
t = naiveBayesClassify(xTest, L$M, L$V, p)
res = evaluate(factor(t[,1]), factor(yTest[,1]))
cat(sprintf('Naive Bayes Classifier\n'))
cat(sprintf('%s (Accuracy)\n', res$accuracy))
cat(sprintf('%s (Precision for class 1)\n', getPrecision(res$confMat, 1)))
cat(sprintf('%s (Recall for class 1)\n', getRecall(res$confMat, 1)))
cat(sprintf('%s (Precision for class 5)\n', getPrecision(res$confMat, 5)))
cat(sprintf('%s (Recall for class 5)\n', getRecall(res$confMat, 5)))
cat('============\n')

### Evaluation Part 3 ###
w0 = c(rep(1,ncol(new_xTrain)))
nIter = 180
w = logisticRegressionWeights(new_xTrain, new_yTrain, w0, nIter)
t = logisticRegressionClassify(new_xTest, w)
res = evaluate(factor(t[,1]), factor(new_yTest[,1]))
cat(sprintf('Logistic Regression Classifier\n'))
cat(sprintf('%s (Accuracy)\n', res$accuracy))
cat(sprintf('%s (Precision for class 0)\n', getPrecision(res$confMat, 0)))
cat(sprintf('%s (Recall for class 0)\n', getRecall(res$confMat, 0)))
cat(sprintf('%s (Precision for class 1)\n', getPrecision(res$confMat, 1)))
cat(sprintf('%s (Recall for class 1)\n', getRecall(res$confMat, 1)))
cat('============\n')

## Evaluation for Naive Bayes on New Ecoli##
p = prior(new_yTrain)
L = likelihood(new_xTrain, new_yTrain)
t = naiveBayesClassify(new_xTest, L$M, L$V, p)
res = evaluate(factor(t[,1]), factor(new_yTest[,1]))
cat(sprintf('Naive Bayes Classifier for New Ecoli\n'))
cat(sprintf('%s (Accuracy)\n', res$accuracy))
cat(sprintf('%s (Precision for class 0)\n', getPrecision(res$confMat, 0)))
cat(sprintf('%s (Recall for class 0)\n', getRecall(res$confMat, 0)))
cat(sprintf('%s (Precision for class 1)\n', getPrecision(res$confMat, 1)))
cat(sprintf('%s (Recall for class 1)\n', getRecall(res$confMat, 1)))
#sink()
