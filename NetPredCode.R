# library needed for the different combination of network analysis and group penalization
library(WGCNA)
library(reshape2)
library(caret)
library(dplyr)
library(tidyr)
library(pROC)
library(grpreg)
library(glasso)
library(MASS)
library(NetComp)
library(mvtnorm)
library(GRridge)
library(SGL)
library(parcor)
library(GeneNet)
library(huge)
library(CompQuadForm)
library(ggm)
library(gglasso)



###############################################################################
#### Functions: Networks, clusters, variable selection ########################
###############################################################################

# Fold creation for cross-validation
# @param y vector containing the ID of each individual in the study.  
# @param k number of folds, default is 10.
# @param list: if True the output is a list containing vectors of indexes for each fold. If Flase,
# the output is a vector of size n (number of individuals) containing the fold index for each individual
# @param nfolds number of folds chosen for the cross-validation.
# @param returnTrain if True the output is a list containing the indexes of the individuals used in the training set for each fold.
createFolds<-function (y, k = 10, list = TRUE, returnTrain = FALSE,seed=seed) 
{
  set.seed(seed)
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2) 
      cuts <- 2
    if (cuts > 5) 
      cuts <- 5
    y <- cut(y, unique(quantile(y, probs = seq(0, 1, length = cuts))), 
             include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      seqVector <- rep(1:k, numInClass[i]%/%k)
      if (numInClass[i]%%k > 0) 
        seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
      foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                        sep = "")
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}

# Double cross-validation for common approaches, lasso, ridge, elastic net
# @param X matrix of features  
# @param Y phenotype values
# @param alpha tuning parameter for elastic net. If alpha =0, ridge regression, and if alpha=1, lasso.
# @param folds list of vector containing the indexes of individuals present in each fold
# @param nfolds number of folds chosen for the cross-validation
# @param family a description of the error distribution and link function to be used in the model.
glmnet.2CV<-function(X,Y,alpha=alpha,folds,nfolds,family=family){
  
  if(nrow(X)!=length(Y)) {stop("Dimensions of X and Y do not match!")}
  alpha=alpha
  n=nrow(X)
  
  
  fit.glmnet=lapply(1:nfolds,function(i)glmnet(X[-folds[[i]],],Y[-folds[[i]]],family=family,standardize=F,alpha=alpha))
  cv.fit.glmnet=lapply(1:nfolds,function(i)cv.glmnet(X[-folds[[i]],],Y[-folds[[i]]],family=family,standardize=F,alpha=alpha))
  
  p.cv.glmnet=unlist(lapply(1:nfolds,function(i)predict(fit.glmnet[[i]],matrix(X[folds[[i]],],ncol=ncol(X)),s=cv.fit.glmnet[[i]]$lambda.min,type="response")))
  p.cv.glmnet<-p.cv.glmnet[order(unlist(folds))]
  
  #Calculate CV mean of the outcome#
  
  
  if (family=="gaussian"){
    p0<-rep(NA,length(Y))
    glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family=family))
    for(i in 1:nfolds){
      p0[folds[[i]]]=rep(coef(glm0[[i]]))}
  }
  else if (family=="binomial"){
    p0<-rep(NA,length(Y))
    glm0<-lapply(1:nfolds,function(i)glm(Y[-folds[[i]]]~1,family=family))
    for(i in 1:nfolds){
      p0[folds[[i]]]= rep(exp(coef(glm0[[i]]))/(1+exp(coef(glm0[[i]]))))}
  }
  
  
  
  
  Q2.glmnet<-(sum((p0-p.cv.glmnet)^2)/sum((Y-p0)^2))
  PRESS.glmnet<-sum((Y-p.cv.glmnet)^2)
  
  return(list(Q2=Q2.glmnet,PRESS=PRESS.glmnet,p=p.cv.glmnet,cv.fit.glmnet=cv.fit.glmnet))
}

# Obtention of the adjacency matrix using graphical lasso
# @param data matrix of features  
GetAdjacencyGLassoHuge=function(data) {
  
  data<-scale(data)
  data_gLasso<-huge(data, lambda = NULL, nlambda = 25, lambda.min.ratio = 0.001, method = "glasso",
                    scr = NULL, scr.num = NULL, cov.output = FALSE, sym = "or", verbose = FALSE)
  model_selection<-huge.select(data_gLasso, criterion = "ric", verbose = FALSE)
  a<-as.matrix(model_selection$opt.icov)
  parameter<-model_selection$opt.lambda
  a[a == "."] <- 0
  a[upper.tri(a)] <- t(a)[upper.tri(a)]
  a<-cov2cor(a)
  adjacency=abs(a)
  results<-list(adjacency,parameter)
  return(results)
}

# Obtention of the adjacency matrix using WGCNA
# @param data matrix of features  
GetAdjacencyWGCNA=function(data){
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  R_sq<-pickSoftThreshold(data, powerVector = powers, verbose = 10,RsquaredCut=0.8)$fitIndices
  power= min(R_sq[which.max(R_sq[,2]),1],R_sq$powerEstimate,min(which(R_sq[,2]>0.8))) 
  adjacency = adjacency(data, power = power,type='unsigned')
  results<-list(adjacency,power)
  return(results)
}

# Clustering analysis based on the dynamic tree cut algorithm
# @param adjacency adjacency matrix of the network studied  
# @param nodes vector of names of the features in the network
GetClusters=function(adjacency,nodes) {
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 10;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  clusters = as.numeric(dynamicMods)
  nodes=data.frame(cbind(nodes,clusters)) # grey - no cluster
  return(nodes)
}

# Clustering analysis based on the dynamic tree cut algorithm
# @param adjacency adjacency matrix of the network studied  
# @param nodes vector of names of the features in the network
FindHubs=function(adjacency,clusters){
  
  hubs=network
  names(hubs)=c("id1","id2","weight")
  names(clusters)=c("nodes","clusters")
  nodes_clustered=clusters[clusters$clusters!=0,]
  hubs=merge(hubs,nodes_clustered,by.x="id1",by.y="nodes")
  names(hubs)[4]="cluster1"
  hubs=merge(hubs,clusters,by.x="id2",by.y="nodes")
  names(hubs)[5]="cluster2"
  
  hubs=hubs %>%
    filter(cluster1==cluster2&weight>0) 
  hubs=data.frame(cbind(c(as.character(hubs$id1),as.character(hubs$id2)),rep(hubs$weight,2),rep(hubs$cluster1,2)))
  names(hubs)=c("id","weight","cluster")
  hubs$weight=as.numeric(hubs$weight)
  
  hubs=hubs %>%
    group_by(id,cluster) %>%
    summarise(
      degree=n(),
      weight=sum(weight)
    )
  hubs=arrange(hubs,cluster,desc(degree),desc(weight))
  
  
  selected_hubs=hubs %>%
    group_by(cluster) %>%
    filter(row_number()<=3) 
  variables=selected_hubs$id    
  
  return(list(variables,hubs,selected_hubs))
}

###############################################################################
################## Functions for Ridge based network-analysis  ################
###############################################################################

EM.mixture <- function(p,eta0,df,tol)
{  # p <- unique partial correlation 
  # eta0, df <- initial values
  # tol <- tolerance
  
  f0 <- function(x,df) {(1-x^2)^((df-3)/2)*(1/beta(1/2,(df-1)/2))}
  fa <- function(x) {dunif(x,-1,1)}
  
  E <- length(p)
  epsilon <- 1000
  i <- 0
  while (epsilon > tol) {
    i <- i+1
    Ak <-  (eta0 * f0(p,df)) / (eta0 * f0(p,df) + (1-eta0) * fa(p)) # conditional expection of the missing 
    expec1 <- (1/2)*(sum(log(1-p^2)*Ak) + (digamma(df/2) - digamma((df-1)/2))*sum(Ak)) # condition expection of first derivative of loglikelihood
    expec2 <- (1/4)*(trigamma(df/2)-trigamma((df-1)/2))*sum(Ak)  # condition expection of second derivative of loglikelihood
    tmp.eta0 <- sum(Ak) / E 
    tmp.df <- df-(1/expec2)*expec1 # using newton-raphson
    epsilon <- max(abs(eta0-tmp.eta0),abs(df-tmp.df))
    eta0 <- tmp.eta0
    df <- tmp.df
  }
  return(list(df=df,eta0=eta0,iter=i))
}

R.separate.ridge <- function(dat,fold,lambda,verbose=FALSE) {
  ## Purpose: Fitting partial correlation matrix using ridge regressions
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dat : p x n data matrix with n (# of sample) and p (# of variables)
  ## -  fold cross validation 
  ## - lambda : the candidate ridge parameter
  ## - verbose : print 
  ## ----------------------------------------------------------------------
  ## Author: Min Jin Ha, Date: 14 September 2013
  
  stopifnot(is.matrix(dat),nrow(dat)>1,ncol(dat)>1)
  n = ncol(dat)
  p = nrow(dat)
  stddat <-(dat - apply(dat,1,mean)) / (apply(dat,1,sd))
  
  coefNI = matrix(0, nrow=p, ncol=p)
  diag.coefNI = lambda.sel = rep(0,p)
  for (i in 1:p){
    if (verbose) cat("variable=", i," ",date(),"\n")
    X       = stddat[-i,]
    y       = stddat[i,]
    fit.lambda = ne.lambda.cv(y=y,XX=X,lambda=lambda,fold=fold) 
    lambdai = fit.lambda$lambda[which.min(fit.lambda$spe)]
    lambda.sel[i] = lambdai
    ridgefit = lm.ridge(y~t(X)-1,lambda=lambdai)
    pred.y = t(X) %*% coef(ridgefit) 
    dvals = svd(X)$d
    diag.coefNI[i] = (n-sum(dvals^2 / (dvals^2 + lambdai)))/sum((y - pred.y)^2)
    coefNI[i,-i] = -diag.coefNI[i]*coef(ridgefit)
  }
  tmp1 <- sqrt(abs(coefNI) * t(abs(coefNI)))
  tmp2 <- sign(coefNI) * upper.tri(coefNI)
  tmp3 <- (tmp2 + t(tmp2)) * tmp1
  diag(tmp3) <- diag.coefNI
  R <- - scaledMat(tmp3)
  diag(R) <- 1
  return(list(R=R,lambda.sel=lambda.sel))
}

Simul.data <- function(G,etaA,N,r,dist.sample="mvnorm")
{
  expression <- array(0,c(G,N,r))
  partialcov <- matrix(rep(0,G*G),ncol=G)
  tri.w <- which(upper.tri(partialcov))
  no.sig <-  ceiling(etaA * length(tri.w))
  sig.node <- sample(tri.w,no.sig)
  partialcov[sig.node] <- runif(no.sig,min=-1,max=1)
  partialcov <- t(partialcov)  + partialcov
  diag(partialcov) <- colSums(abs(partialcov))+0.0001  # to make p.d. matrix
  partialcor <- partialcov / sqrt(diag(partialcov) %*% t(diag(partialcov)))
  invcor <- - partialcor
  diag(invcor) <- 1
  truecor <- solve(invcor)
  truecor.scaled=truecor/ sqrt(diag(truecor) %*% t(diag(truecor)))
  if (dist.sample == "mvnorm") {
    for (i in 1:r)  {
      expression[,,i] <- t(rmvnorm(N,mean=rep(0,G),truecor))
    }
  }else if (dist.sample == "mvt") {
    for (i in 1:r)  {
      expression[,,i] <- t(rmvt(n=N,sigma=truecor,df=2))
    }
  }
  return(list(simul.data=expression,true.partialcor=partialcor,
              truecor.scaled=truecor.scaled,sig.node=sig.node))
}

Structured.estimate <- function(X,E) {
  # Using rank min(n,p) inverse
  # X <- std data matrix
  # E estimated edge set
  E <- matrix(E,ncol=2)
  p <- nrow(X)
  n <- ncol(X)
  index <- 1:p
  K <- matrix(0,p,p)
  RSS <- rep(0,p)
  diag(K) <- apply(X,1,var)
  for (response.id in 1:p) {
    y <- X[response.id,] 
    ne <- c(E[which(E[,1]==response.id),2],E[which(E[,2]==response.id),1])
    if (length(ne) == 0) RSS[response.id] <- sum(y^2)
    if (length(ne) > 0){
      D <- t(matrix(X[ne,],ncol=n)) 
      svdD = svd(D)
      coef <- svdD$v%*%(t(svdD$v)/svdD$d^2) %*% t(D) %*% y
      RSS[response.id] <- sum((y - D %*% coef)^2)
      resid.var <- RSS[response.id] /(n-length(ne)) 
      coef.cl <- c(1/resid.var,-coef/resid.var)
      id.cl <- c(response.id,ne)  
      o.cl <- order(id.cl)
      K[response.id,id.cl[o.cl]] <- coef.cl[o.cl]
    }
  }
  tmp1 <- sqrt(abs(K) * t(abs(K)))
  tmp2 <- sign(K) * upper.tri(K)
  tmp3 <- (tmp2 + t(tmp2)) * tmp1
  diag(tmp3) <- diag(K)
  R <- - scaledMat(tmp3)
  diag(R) <- 1
  return(list(R=R,K=tmp3,RSS=RSS))
}

getEfronp <- function(z,num.bin=120,max.q=9,pct=0,pct0=0.25,cc=1.2,plotIt=F) {
  # num.bin <- number of bins, max.q <- maximum degree of polynomial considered 
  # pct :   Proportion of the z distribution used in fitting the marginal density f(z) by polynomial regression.
  
  if (num.bin<2) stop("the number of bins has to be greater than 1") 
  if (pct>pct0) stop("pch0 must be greater than or equal to pct")
  stopifnot(pct>=0)
  
  pvalues <- 2*(1-pnorm(abs(z)))
  
  # Z values for marginal density estimation
  v = quantile(z, c(pct, 1 - pct))
  lo = v[1]
  hi = v[2]
  zz = pmax(pmin(z, hi), lo)
  
  # X values for marginal density estimation
  breaks = seq(lo,hi,length=num.bin+1)
  hcomb = hist(zz,breaks=breaks,plot=FALSE) 
  x = hcomb$mids
  s = hcomb$counts
  X <- sapply(1:max.q, function(j) x^j)
  
  kstest = mu0hat = sigma0hat = eta = rep(NA,max.q)
  kstest[1] = ks.test(pvalues,"punif")$statistic
  mu0hat[1] = 0
  sigma0hat[1] = 1
  
  cof <- 0
  for (p in 2:max.q){
    tryCatch({cof <- coef(glm(s ~ X[,1:p],family="poisson"))},error = function(e) cof=0)
    if (length(cof)==p+1 & sum(is.na(cof))==0 ){
      fz <- exp(sapply(0:p,function(i) zz^i)%*%cof)
      mu0hat[p] <- z[which.max(fz)]
      tmp <- sapply(0:(p-2),function(i)(i+2)*(i+1)*mu0hat[p]^i)
      logfzdiff2 <- t(as.matrix(tmp)) %*% as.matrix(cof[-(1:2)])
      if (logfzdiff2 < 0) {
        sigma0hat[p] <- (-logfzdiff2)^(-1/2)
        correctZ <- (z-mu0hat[p])/sigma0hat[p]
        correctp <- 2*(1-pnorm(abs(correctZ)))
        kstest[p] <- ks.stat.p(correctp)
      }  
    }
  }
  
  qq <- which.min(kstest[-1])+1
  correctZ <- (z-mu0hat[qq])/sigma0hat[qq]
  correctp <- 2*(1-pnorm(abs(correctZ)))
  mu = mu0hat[qq]
  sigma = sigma0hat[qq]
  
  # Eta calculation (null proportion)
  fit = glm(s ~ X[,1:qq], poisson)
  f = fit$fit
  l = log(f)
  xmax = x[which.max(l)]
  
  v0 = quantile(z, c(pct0, 1 - pct0))
  lo0 = v0[1]
  hi0 = v0[2]
  idx0 = which(x>lo0 & x<hi0)
  # cat(length(idx0),"\n")
  if (length(idx0) > 10) {
    x0 = x[idx0]
    y0 = l[idx0]
    X0 <- cbind(x0 - mu, (x0 - mu)^2)
    fit.null <- lm(y0 ~ X0)
    coef.null <- fit.null$coef
    #X00 <- cbind(1, x - mu, (x - mu)^2)
    #l0 <- as.vector(X00 %*% coef.null)
    #f0 <- exp(l0)
    #eta <- max(0,min(sum(f0)/sum(f),1))
  }else {
    #cat("too small number of x for null region\n")
    idz0 = which(z>lo0 & z<hi0)
    z0 = z[idz0]
    fz = exp(sapply(0:qq,function(i) z^i)%*%coef(fit))
    lz = log(fz)
    y0 = lz[idz0]
    Z0 = cbind(z0-mu,(z0-mu)^2)
    fit.null = lm(y0 ~ Z0)
    coef.null = fit.null$coef
    # Z00 = cbind(1,x-mu,(x-mu)^2) 
    # l0 = as.vector(Z00 %*% coef.null)
    # f0 = exp(l0)
    # eta = max(0,min(sum(f0)/sum(f),1))
  }
  ### Correct mu and sigma
  mu = -coef.null[2]/(2*coef.null[3]) + mu
  sigma = sqrt(-1/(2*coef.null[3]))
  lo00 = mu - cc*sigma
  hi00 = mu + cc*sigma
  pi = sum(z>lo00 & z<hi00) / length(z)
  G0 = 1- 2*pnorm(lo00,mean=mu,sd=sigma)
  eta = max(min(pi / G0,1),0)
  
  if (plotIt) {
    oz = order(zz)
    h = hist(zz,breaks=breaks,xlab=expression(psi),main="")
    lines(x,f,col=3,lwd=2)
    lines(zz[oz],dnorm(zz,mu,sigma)[oz]*length(zz)*mean(diff(h$breaks)),col=2,lty=2,lwd=1.5)
    legend("topright",legend=c("f",expression(f[0])),col=c(3,2),lty=1:2,lwd=c(2,1.5),bty="n",seg.len=3)
  }
  
  return(list(correctz=correctZ,correctp=correctp,q=qq,mu0hat=mu,sigma0hat=sigma,etahat=eta))
}

ks.stat.p<-function(p){
  ### This function calculates kolmogorov-smirnov statistic for p-values 
  ### Compare p values with uniform distribution
  ### Modified from ks.test function
  ### Input
  ### - p : a numeric vector indicating pvalues
  y <- get("punif", mode = "function", envir = parent.frame())
  p <- p[!is.na(p)]
  n <- length(p)
  TIES <- FALSE
  if (length(unique(p)) < n) TIES <- TRUE
  p <- y(sort(p)) - (0:(n - 1))/n
  STATISTIC <- max(c(p, 1/n - p))
  STATISTIC
}

lambda.TargetD<-function(X) {
  # X <- centered data for covariance shrinkage and standardized data for correlation shrinkage
  G <- nrow(X); N <-ncol(X)
  wmean <- X %*% t(X)/N
  S <- (N/(N-1)) * wmean  
  varw <- (( X^2 %*% t(X^2) ) - (N*wmean^2)) * (N/((N-1)^3))
  diag.w <- which(row(wmean)==col(wmean))
  lambda <- max(0,min(sum(varw[-diag.w]) / sum((S^2)[-diag.w]),1))
  return(lambda)
}

lambda.cv<-function(X,lambda,fold) {
  ## X <- p by n
  ## fold-cross validation
  ## lambda <- XX scale
  p <- nrow(X); n <- ncol(X)
  cv <- 1:n%%fold+1
  cv <- cbind(cv,sample(1:n))
  spe <- rep(0,length(lambda))
  for (cv.index in 1:fold) {
    testindex <- cv[cv[,1]==cv.index,2]
    train <- matrix(X[,-testindex],nrow=p)
    test <- matrix(X[,testindex],nrow=p)
    std.train <- (train - apply(train,1,mean)) / (apply(train,1,sd))
    if (ncol(test)==1) std.test <- test
    if (ncol(test)>1) std.test <- (test - apply(test,1,mean)) / (apply(test,1,sd))
    k <- length(lambda) # no. of candidate lambda
    spe.mat <- matrix(nrow=p,ncol=k)
    for (response.id in 1:p) {
      y <- std.train[response.id,] # train response
      D <- t(std.train[-response.id,]) # train desing matrix (sample by variables)
      #coef1 <- coef(lm.ridge(y~D-1,lambda=lambda.array)) # lm.ridge used
      Ds <- svd(D) # svd of D
      d <- Ds$d # singular values of D
      dx <- length(d) # min(n,p)
      div <- d^2 + rep(lambda, rep(dx, k))
      rhs <- t(Ds$u) %*% y
      a <- drop(d * rhs)/div
      dim(a) <- c(dx, k)
      coef <- Ds$v %*% a
      y.test <- std.test[response.id,]
      D.test <- t(std.test[-response.id,])
      spe.mat[response.id,] <- colSums((y.test - D.test %*% coef)^2) 
    }
    spe <- spe + colSums(spe.mat)
  }
  return(list(lambda=lambda,spe=spe/(fold*p)))
}

lambda.pcut.cv <- function(X,lambda,pcut,fold=10) {
  ## X <- p by n
  ## fold-cross validation
  ## lambda <- XX scale
  p <- nrow(X); n <- ncol(X)
  cv <- 1:n%%fold+1
  cv <- cbind(cv,sample(1:n))
  
  PE = matrix(0,nrow=length(lambda),ncol=length(pcut))
  rownames(PE) = lambda
  colnames(PE) = pcut
  cat("cv:")
  for (cv.index in 1:fold) {
    cat(" ", cv.index)
    testindex <- cv[cv[,1]==cv.index,2]
    train <- matrix(X[,-testindex],nrow=p)
    test <- matrix(X[,testindex],nrow=p)
    std.train <- (train - apply(train, 1, mean))/(apply(train, 1, sd))
    
    if (ncol(test)==1) std.test <- test
    if (ncol(test)>1) std.test <- (test - apply(test,1,mean)) / (apply(test,1,sd))
    k <- length(lambda) # no. of candidate lambda
    spe.mat <- matrix(nrow=p,ncol=k)
    
    PE = PE + lambda.pcut.cv1(train=std.train,test=std.test,lambda,pcut)
  }
  return(PE/fold)
}

lambda.pcut.cv1 <- function(train,test,lambda,pcut){
  # X <- std data
  p <- nrow(test)
  k <- length(lambda)
  w.upper = which(upper.tri(diag(p)))
  w.array = which(upper.tri(diag(p)),arr.ind=T)
  tunings <- matrix(ncol=2,nrow=length(lambda) * length(pcut))
  S <- cor(t(train))
  R <- matrix(nrow=length(w.upper),ncol=k)
  #diag.W <- matrix(nrow=p,ncol=k)
  kk <- 0
  for (la in lambda) {
    kk <- kk+1
    tmp <- solve(S + la * diag(p))
    # diag.W[,kk] <- diag(tmp)
    R[,kk] <- -scaledMat(tmp)[w.upper]
  }
  transR <- trans.Fisher(R)
  efron <- sapply(1:k,function(i) getEfronp(transR[,i])$correctp)
  
  risk <- matrix(nrow=k,ncol=length(pcut))
  colnames(risk) <- pcut
  rownames(risk) <- lambda
  for (th in pcut) {
    thR <- R * (efron < th)
    coef.mat <- matrix(nrow=p*p,ncol=k)
    for (kk in 1:k) {
      E = w.array[thR[,kk]!=0]
      fit <- Structured.estimate(train,E)     
      coef.mat[,kk] <- c(diag(sqrt(diag(fit$K))) %*% fit$R %*% diag(1/sqrt(diag(fit$K)))) # p by p
    }
    
    CV <- 0
    for (response.id in 1:p) {
      n = ncol(test)
      y <- test[response.id,] 
      D <- t(test[-response.id,]) 
      #I_P <- sapply(lambda,function(la) diag(n) - (D %*% solve(t(D) %*% D + la*diag(p-1)) %*% t(D)))
      #Pinv <- sapply(1:length(lambda),function(i) diag(1/diag(matrix(I_P[,i],ncol=n))))
      #loss <- sapply(1:length(lambda), function(i) (matrix(Pinv[,i],ncol=n)%*%(y - D %*% (coef.mat[((response.id-1)*p+1):(response.id*p) ,i])[-1]))^2)
      loss <- sapply(1:length(lambda), function(i)(y - D %*% (coef.mat[((response.id-1)*p+1):(response.id*p) ,i])[-1])^2)
      CV <- CV + colSums(loss) 
    }
    risk[,which(pcut==th)] <- CV / (n*p) 
  }
  #w <- which(risk==min(risk), arr.ind=T)
  #setlambda <- lambda[w[1,1]] 
  #setth <- pcut[w[1,2]]
  #thR <- R[,which(lambda==setlambda)] * (efron[,which(lambda==setlambda)] < setth)
  #mat.thR <- matrix(0,nrow=p,ncol=p)
  #mat.thR[w.upper] <- thR
  #mat.thR <- mat.thR + t(mat.thR) 
  #return(list(lambda=lambda,pcut=pcut,risk=risk,thR = mat.thR,setlambda=setlambda,setth=setth))
  return(risk)
}

ne.lambda.cv <- function(y,XX,lambda,fold) {
  # y : p vector
  # XX: p by n 
  # fold cross validation
  
  X = rbind(y,XX)
  response.id = 1
  p = nrow(X); n <- ncol(X)
  cv = 1:n%%fold+1
  cv = cbind(cv,sample(1:n))
  k = length(lambda)
  spe = rep(0,k)
  for (cv.index in 1:fold) {
    testindex <- cv[cv[,1]==cv.index,2]
    train <- matrix(X[,-testindex],nrow=p)
    test <- matrix(X[,testindex],nrow=p)
    std.train <- (train - apply(train,1,mean)) / (apply(train,1,sd))
    if (ncol(test)==1) std.test <- test
    if (ncol(test)>1) std.test <- (test - apply(test,1,mean)) / (apply(test,1,sd))
    
    y <- std.train[response.id,] # train response
    D <- t(std.train[-response.id,]) # train desing matrix (sample by variables)
    #coef1 <- coef(lm.ridge(y~D-1,lambda=lambda.array)) # lm.ridge used
    Ds <- svd(D) # svd of D
    d <- Ds$d # singular values of D
    dx <- length(d) # min(n,p)
    div <- d^2 + rep(lambda, rep(dx, k))
    rhs <- t(Ds$u) %*% y
    a <- drop(d * rhs)/div
    dim(a) <- c(dx, k)
    coef <- Ds$v %*% a
    y.test <- std.test[response.id,]
    D.test <- t(std.test[-response.id,])
    spe <- spe + colSums((y.test - D.test %*% coef)^2) 
  }
  return(list(lambda=lambda,spe=spe/fold))
}

scaledMat <- function(x){
  newx=x/sqrt(diag(x) %*% t(diag(x)))
  return(newx)
}

trans.Fisher <-  function(x) {
  x[x>=(1-1e-07)] <- 1 - 1e-07
  x[x<=(-1+1e-07)] <- -1 + 1e-07
  return(log((1+x)/(1-x))/2)
}

datasetForVisualization<-function(sparseMatrix,names){
  data<-NULL
  for (i in 1:(length(sparseMatrix[1,])-1)){
    nodeConnections<-cbind(rep(names[i],length(sparseMatrix[1,])-i),names[c((i+1):length(names))],sparseMatrix[i,c((i+1):length(names))])
    data<-rbind(data,nodeConnections)
  }
  rownames(data)<-c(1:length(data[,1]))
  data
}

#Standardization of each glycan
standardization<-function(matrix){
  data<-NULL
  for (i in 1:length(names(matrix))){
    standard<-(matrix[,i]-mean(matrix[,i]))/sd(matrix[,i])
    data<-cbind(data,standard)
  }
  return(data)
}

#Function to obtain networks based on partial correlation using ridge penalty
#parameters dataset and threshold for partial correlation testing
GetAdjacencyHa<-function(data,th){
  
  ##########################Parameters needed for applying Ha's approach
  n<-length(data[,1])
  p<-length(data[1,])
  w.upper = which(upper.tri(diag(p)))
  w.array = which(upper.tri(diag(p)),arr.ind=TRUE)
  
  ################################ Ha's approach 
  data<-scale(data,center=T,scale=T)
  
  # estimate ridge parameter
  lambda.array = seq(0.005,20,by=0.1) * (n-1)
  fit = lambda.cv(t(data),lambda.array,fold=10)
  lambda = fit$lambda[which.min(fit$spe)]/(n-1)
  
  # calculate partial correlation using ridge inverse
  partial = solve(lambda*diag(p) + cor((data)))
  partial = (-scaledMat(partial))[w.upper]
  
  # get p-values from empirical null distribution 
  efron.fit <- getEfronp(trans.Fisher(partial),num.bin=50,max.q=13)
  
  # estimate the edge set of partial correlation graph with FDR control at level 0.01
  wsig <- which(p.adjust(efron.fit$correctp,method="bonferroni") < th )
  E <- w.array[wsig,]
  
  # structured estimation     
  fit <- Structured.estimate(t(data),E)
  th.partial <- fit$R
  adjacency<-abs(th.partial)
  results<-list(adjacency,lambda)
  return(results)
}


###############################################################################
################## Main function  #############################################
###############################################################################

# Double cross-validation function for the pred net 
# @param method1 method used for network analysis possible values are "WGCNA" "gLasso" or "Ridge"
# @param method2 method used for prediction model possible values are "hubs" "SGL" or "grpLasso" "gRidge" 
# @param data matrix or dataframe of features to built the prediction model
# @param phenotype vector containing outcome's values
# @param nfoldsOut number of fold for the cross-validation outer loop
# @param nfoldsIn number of fold the cross-validation inner loop
# @param alpha parameter for the SGL method only
# @param seed 
NetPred2CV<-function(method1,method2,data,phenotype,nfoldsOut,nfoldsIn,alpha,seed){
  X<-data
  outcome<-phenotype
  names(data)<-c(1:length(data))
  options(stringsAsFactors = FALSE)
  N=dim(data)[2]
  nodes=c(1:ncol(data))
  
  folds=createFolds(1:nrow(data), k = nfoldsOut, list = T,returnTrain=F,seed=seed)
  coef_cv=NULL
  predict.fit=NULL
  modules=NULL
  ncluster=NULL
  for (i in 1:nfoldsOut){
    data1=data[-folds[[i]],] # train data
    phenotype1=phenotype[-folds[[i]]]
    data.val<-data[folds[[i]],] # test data
    
    #select method of network construction
    if (method1=="WGCNA"){
      adjacency=GetAdjacencyWGCNA(data1)[[1]]
    } else if (method1=="gLasso") {
      adjacency=GetAdjacencyGLassoHuge(data1)[[1]]
    } else if (method1=="Ridge"){
      adjacency=GetAdjacencyHa(data1,0.05)[[1]]
    }
    #apply clusterisation algorythm
    clusters=GetClusters(adjacency,nodes)
    modules<-cbind(modules,clusters[,2])
    ncluster<-c(ncluster,length(unique(clusters[,2])))
    clusters0=as.numeric(clusters[,2])
    clusters0[which(clusters0==0)]=max(clusters0)+1
    network=adjacency
    network[upper.tri(network)]=NA
    diag(network)=NA
    network=melt(network,na.rm = TRUE)
    
    #choose method of variable selection
    if (method2=="hubs"){
      variablesHubs<-rep(0,length(X[1,])) 
      if (length(unique(clusters[,2]))==1){
        fit.ridge<-glmnet(data.matrix(X[-folds[[i]],]),phenotype1,family="gaussian",standardize=F,alpha=0)
        cv.fit.ridge<-cv.glmnet(data.matrix(X[-folds[[i]],]),phenotype1,family="gaussian",standardize=F,alpha=0,nfolds=nfoldsIn)
        predict.fit[[i]]<-predict(fit.ridge,data.matrix(data.val),s=cv.fit.ridge$lambda.min,type="response")
        coef_cv[[i]]<-coef(fit.ridge,as.matrix(data.val),s=cv.fit.ridge$lambda.min)[-1]
      }
      else{hubs=FindHubs(network,clusters)[[1]] # one hub per cluster
      hubs<-as.numeric(hubs)
      variables=X[-folds[[i]],hubs] # subset of data to use in prediction
      variables.val=data.val[,hubs]
      variables_cv<-hubs # list of variables selected on each iteration
      fit.glm <- glmnet(data.matrix(variables),phenotype1,family="gaussian",standardize=F,alpha=0)
      cv.fit.glm <- cv.glmnet(data.matrix(variables),phenotype1,family="gaussian",standardize=F,alpha=0,nfolds=nfoldsIn)
      coef_cvHubs<-coef(fit.glm ,as.matrix(variables.val),s=cv.fit.glm$lambda.min)[-1]
      predict.fit[[i]] <- predict(fit.glm,data.matrix(variables.val),s=cv.fit.glm$lambda.min,type="response")
      for (o in 1:length(hubs)){
        index<-hubs[o]
        variablesHubs[index]<-coef_cvHubs[o]
      }
      coef_cv[[i]]<-variablesHubs
      }
    }
    
    if (method2=="SGL"){
      dataSGL<-list(x=data1,y=phenotype1)
      cv.fit.sglm<-cvSGL(dataSGL,clusters0,type="linear",nlam=20,nfold=nfoldsIn,alpha=alpha)
      fit.sglm <-SGL(dataSGL,clusters0,type="linear",nlam=50,alpha=alpha)
      predict.fit[[i]] <- predictSGL(x=fit.sglm,newX=as.matrix(data.val),lam=which.min(abs(fit.sglm$lambdas-cv.fit.sglm$lambdas[which.min(cv.fit.sglm$lldiff)])))
      coef_cv[[i]]<-fit.sglm$beta[,which.min(abs(fit.sglm$lambdas-cv.fit.sglm$lambdas[which.min(cv.fit.sglm$lldiff)]))]
    }
    
    
    if (method2=="grplasso"){    
      fit.gglasso=gglasso(as.matrix(data1),phenotype1,nlam=50,lambda.factor=0.01,loss="ls",group=clusters0)
      cv.fit.gglasso=cv.gglasso(as.matrix(data1),phenotype1,lambda.factor=0.01,nlam=50,loss="ls",group=clusters0,nfolds=,nfolds=nfoldsIn)
      predict.fit[[i]]=unlist(predict(fit.gglasso,as.matrix(data.val),s=cv.fit.gglasso$lambda.min,type="link"))
      coef_cv[[i]]<- coef(fit.gglasso,as.matrix(data.val),s=cv.fit.gglasso$lambda.min)[-1]
    }
    
    if (method2=="gRidge"){     
      clusters1<-as.factor(clusters0)
      partition_ridge<-CreatePartition(clusters1)
      fit.grridge=try(grridge(t(data.matrix(data1)),phenotype1,partitions=list(partition_ridge)))
      predict.fit[[i]]=try(predict.grridge(fit.grridge,t(data.matrix(data.val)))[,2])
      coef_cv[[i]]<-try(fit.grridge$betas)
    }
  }
  ###################################################
  p0<-rep(NA,length(phenotype))
  glm0<-lapply(1:length(folds),function(i)glm(phenotype[-folds[[i]]]~1,family="gaussian"))
  for(i in 1:length(folds)){
    p0[folds[[i]]]=rep(coef(glm0[[i]]))
  }
  
  predict.fit.unlisted=unlist(predict.fit)
  predict.fit.unlisted <- predict.fit.unlisted[order(unlist(folds))]
  Q2.new<-(sum((predict.fit.unlisted-p0)^2)/sum((phenotype-p0)^2))
  PRESS<-sum((predict.fit.unlisted-phenotype)^2)
  
  Ridge<-glmnet.2CV(X=data.matrix(X),Y=outcome,alpha=0,folds=folds,nfolds=nfolds,family="gaussian");
  Q2.Ridge<-Ridge$Q2;
  PRESSRIDGE<-Ridge$PRESS
  
  Lasso<-glmnet.2CV(X=data.matrix(X),Y=outcome,alpha=1,folds=folds,nfolds=nfolds,family="gaussian");
  Q2.Lasso<-Lasso$Q2;
  PRESSLASSO<-Lasso$PRESS
  
  ElasticNet<-glmnet.2CV(X=data.matrix(X),Y=outcome,alpha=0.5,folds=folds,nfolds=nfolds,family="gaussian");
  Q2.ElasticNet<-ElasticNet$Q2;
  PRESSElasticNet<-ElasticNet$PRESS
  
  try(cat("Q2 & PRESS:" ,Q2.new,PRESS," ,
          Results Lasso:",Q2.Lasso,PRESSLASSO," ,
          Results Elastic Net:",Q2.ElasticNet,PRESSElasticNet,",
          Results Ridge:",Q2.Ridge,  PRESSRIDGE))
  
  coefficientcv<-Reduce(cbind,coef_cv)
  results<-list(predAc=c(Q2.new,PRESS),coef=coefficientcv,predict.fit.unlisted)
  results
}

