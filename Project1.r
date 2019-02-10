##################################
###### STAT 557 (Project 1) ######
##################################

rm(list=ls()) ## To clear your environment

## Read the data
xTrain=read.csv("ecoli_xTrain.csv")
yTrain=read.csv("ecoli_yTrain.csv")
xTest=read.csv("ecoli_xTest.csv")
yTest=read.csv("ecoli_yTest.csv")


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
    df = cbind(xTrain, class=yTrain$X1)
    M = matrix(0, f, c)
    V = matrix(0, f, c)
    for (j in 1:c) {
        df1 = df[df$class==j, 1:f] # Get the observations for class-j
        M[,j] = as.numeric(apply(df1, 2, mean)) # Get mean for each feature
        V[,j] = as.numeric(apply(df1, 2, var)) # Get var for each feature
    }
    return(list(M=M, V=V))
}

naiveBayesClassify <- function(xTest, M, V, p){
    feature = nrow(M)
    cls = ncol(M)
    rows = nrow(xTest)
    t = cbind(c(rep(0, rows)))
    for (i in 1:rows) {
        pst = c(rep(0,cls))
        for (c in 1:cls) {
            llh = c(rep(0,feature))
            for (f in 1:feature) {
                llh[f] = dnorm(xTest[i,f], M[f,c], sqrt(V[f,c]), log)
            }
            pst[c] = logProd(c(llh, log(p[c])))
        }
        t[i] = which.max(pst)
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
    return(t)
}
