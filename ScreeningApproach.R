##########################################
##
##  Reproduce Wanghuan's Method for feature screening
##  Written by Junli Lin
##
##########################################

# isCluster = TRUE # in cluster
if (exists("isCluster"))
{
  if (isCluster)
  {
    lib.loc = "~/work/R/libraries"
  }
}

if (!exists("lib.loc"))
{
  lib.loc = .libPaths()
}

library(MASS) # for mvrnorm()
library(splines) # for bs()
library(parallel)
library(iterators, lib.loc = lib.loc) # for gee()
library(gee, lib.loc = lib.loc) # for gee()
library(expm, lib.loc = lib.loc) # for sqrtm()
library(foreach, lib.loc = lib.loc)
library(snow, lib.loc = lib.loc)
library(doSNOW, lib.loc = lib.loc)

BETA.CONSTANT = c(0.4, 0.5, -0.4, -0.5)

## Generate all tij that will be used
genTime = function(N){
  tmp = rnorm(N)
  sort(pnorm((tmp - mean(tmp))/sd(tmp)))
}

## Generate the covariance structure
genCovMat = function(t.vec, rho1 = 0.6, rho2 = 0.4){
  n = length(t.vec)
  var.vec = 0.5 + 3 * t.vec^3
  mat1 = matrix(1:n, n, n)
  cor.mat = 0.5 * rho1^abs(mat1 - t(mat1)) + 0.5 * rho2
  diag(cor.mat) = 1
  cov.mat = diag(sqrt(var.vec)) %*% cor.mat %*% diag(sqrt(var.vec))
  cov.mat
}

## Generate epsilons from the above covariance structure.
genEpsilon = function(t.vec) {
  n = length(t.vec)
  mvrnorm(n = 1, mu = rep(0, n), Sigma = genCovMat(t.vec))
}

## Generate values of the gamma functions at all times
## In this case, they are fixed functions.
## p is the number of SNP
genGamma = function(t.vec, p) {
  n = length(t.vec)
  significant.vec = BETA.CONSTANT
  matrix(c(significant.vec, rep(0, p-4)), n, p, byrow = T)
}

## Generate values of the beta functions at all times
## In this case, they are fixed functions
## q is t he number of non0intercept covariates
genBeta = function(t.vec, q) {
  n = length(t.vec)
  matrix((0:q+1)/10, n, q+1, byrow=T)
}

## Generate values of the gender
## In this case, they are randomly 0 or 1.
genGender = function(n.sub) {
  sample(0:1, n.sub, replace = T)
}

## Generate values of the SNP
## In this case, they are randomly 0, 1 or 2.
genSNP = function(n.sub, p) {
  matrix(sample(0:2, n.sub*p, replace = T), n.sub, p)
}

## Generate data
genSimulatedData = function(n.sub, n.rep, p)
{
  id = 1:n.sub
  id.all = rep(1:n.sub, each = n.rep)
  gender = genGender(n.sub)
  gender.all = gender[id.all]
  SNP = genSNP(n.sub, p)
  SNP.all = SNP[id.all, ]
  
  time.all = genTime(n.sub * n.rep)
  beta.all = genBeta(time.all, 1)
  gamma.all = genGamma(time.all, p)
  epsilon.all = rep(NA, n.sub * n.rep)
  for (i in 1:n.sub) {
    index.current = (id.all == i)
    epsilon.all[index.current] = genEpsilon(time.all[index.current])
  }
  y.all = rowSums(beta.all * cbind(1, gender.all)) + rowSums(gamma.all * SNP.all) + epsilon.all
  
  simData = data.frame(id.all, y.all, time.all, gender.all, SNP.all)
  simData
}

# Step 1: Construct the working model with only the kth component
# Step 2: Construct the variance structure of the working model (with kth component)
# Step 3: Use the variance matrix and WLS to compute the coefficients for the working model
# Step 4: Compute the residuals and the screening criterion

wanghuanMethod = function(y.all, time.all, gender.all, SNP.all, n.sub, n.rep, 
                          degree.spline = 3, knot.spline = 4, Mv = 15) {
  
  # basic dimensions
  p = ncol(SNP.all)
  id.all = rep(1:n.sub, each = n.rep)
  
  # generate b-spline functions for all the times
  total.df.spline = 1 + knot.spline + degree.spline
  time.spline = cbind(1, bs(time.all, df = total.df.spline - 1, degree = degree.spline, intercept = F))
  
  # generate the design matrix for non-SNP predictors
  nonSNP.repeated = cbind(1, gender.all)[, rep(1:2, each = total.df.spline)]
  time.spline.repeated = time.spline[, rep(1:total.df.spline, 2)]
  nonSNP.design.mat = nonSNP.repeated * time.spline.repeated
  
  # generate the design matrix for SNP predictors
  SNP.repeated = SNP.all[, rep(1:p, each = total.df.spline)]
  time.spline.repeated2 = time.spline[, rep(1:total.df.spline, p)]
  SNP.design.mat = SNP.repeated * time.spline.repeated2
  
  # fit model (6) in the paper assuming i.i.d. error
  model6 = lm(y.all ~ 0 + nonSNP.design.mat)
  var.estimated.all = fitted(lm(model6$residuals^2 ~ 0 + time.spline)) # MAY BE NEGATIVE??
  
  # divide the estimated variances from model (6) so that the variances are all 1
  y.all.1 = y.all / sqrt(var.estimated.all)
  nonSNP.design.mat.1 = nonSNP.design.mat / sqrt(var.estimated.all)
  
  # fit model (6) in the paper assuming stationary M-dependent correlation structures,
  # the gee function requires that the numbers of observations in the groups are the same
  # check this!
  model6.2 = gee(y.all.1 ~ 0 + nonSNP.design.mat.1, id = id.all, family = gaussian("identity"), corstr ="stat_M_dep", Mv = Mv)
  cor.mat.estimated = model6.2$working.correlation
  cor.mat.estimated.inv = solve(cor.mat.estimated)
  
  # store the individual weight matrixes
  W.0.5.all = matrix(NA, n.rep*n.sub, n.rep)
  for (i in 1:n.sub) {
    id.i = (id.all == i)
    var.estimated.i = var.estimated.all[id.i]
    var.estimated.i.inv.5.mat = diag(1/sqrt(var.estimated.i))
    W.i = var.estimated.i.inv.5.mat %*% cor.mat.estimated.inv %*% var.estimated.i.inv.5.mat / n.rep
    W.0.5.all[id.i, ] = sqrtm(W.i)
  }
  
  # fit model (5) in the paper for each kth SNP
  u.vec = rep(NA, p)
  for (k in 1:p) {
    pred = as.matrix(cbind(nonSNP.design.mat, SNP.design.mat[(k-1)*total.df.spline+1:total.df.spline]))
    resp = y.all
    for (i in 1:n.sub) {
      id.i = (id.all == i)
      W.0.5.i = W.0.5.all[(i-1)*n.rep + 1:n.rep, ]
      pred[id.i, ] = W.0.5.i %*% pred[id.i, ]
      resp[id.i] = W.0.5.i %*% resp[id.i]
    }
    resp.fitted = fitted(lm(resp ~ 0 + pred))
    u.vec[k] = mean((resp.fitted - resp)^2)
  }
  
  u.vec
}



## For this version, SNP.all does not have the same number of rows as others.
## Only 1 row in SNP.all is associated with each subject.

wanghuanMethod.2 = function(y.all, time.all, gender.all, SNP.all.2, n.sub, n.rep, 
                          degree.spline = 3, knot.spline = 4, Mv = 15) {

  # basic dimensions
  p = ncol(SNP.all.2)
  id.all = rep(1:n.sub, each = n.rep)
  
  # generate b-spline functions for all the times
  total.df.spline = 1 + knot.spline + degree.spline
  time.spline = cbind(1, bs(time.all, df = total.df.spline - 1, degree = degree.spline, intercept = F))
  
  # generate the design matrix for non-SNP predictors
  nonSNP.repeated = cbind(1, gender.all)[, rep(1:2, each = total.df.spline)]
  time.spline.repeated = time.spline[, rep(1:total.df.spline, 2)]
  nonSNP.design.mat = nonSNP.repeated * time.spline.repeated
  
  # generate the design matrix for SNP predictors
  SNP.repeated = SNP.all.2[rep(1:n.sub, each = n.rep), rep(1:p, each = total.df.spline)]
  time.spline.repeated2 = time.spline[, rep(1:total.df.spline, p)]
  SNP.design.mat = SNP.repeated * time.spline.repeated2
  
  # fit model (6) in the paper assuming i.i.d. error
  model6 = lm(y.all ~ 0 + nonSNP.design.mat)
  var.estimated.all = fitted(lm(model6$residuals^2 ~ 0 + time.spline)) # MAY BE NEGATIVE??
  
  # divide the estimated variances from model (6) so that the variances are all 1
  y.all.1 = y.all / sqrt(var.estimated.all)
  nonSNP.design.mat.1 = nonSNP.design.mat / sqrt(var.estimated.all)
  
  # fit model (6) in the paper assuming stationary M-dependent correlation structures,
  # the gee function requires that the numbers of observations in the groups are the same
  # check this!
  model6.2 = gee(y.all.1 ~ 0 + nonSNP.design.mat.1, id = id.all, family = gaussian("identity"), corstr ="stat_M_dep", Mv = Mv)
  cor.mat.estimated = model6.2$working.correlation
  cor.mat.estimated.inv = solve(cor.mat.estimated)
  
  # store the individual weight matrixes
  W.0.5.all = matrix(NA, n.rep*n.sub, n.rep)
  for (i in 1:n.sub) {
    id.i = (id.all == i)
    var.estimated.i = var.estimated.all[id.i]
    var.estimated.i.inv.5.mat = diag(1/sqrt(var.estimated.i))
    W.i = var.estimated.i.inv.5.mat %*% cor.mat.estimated.inv %*% var.estimated.i.inv.5.mat / n.rep
    W.0.5.all[id.i, ] = sqrtm(W.i)
  }
  
  # fit model (5) in the paper for each kth SNP
  u.vec = rep(NA, p)
  for (k in 1:p) {
    current.SNP.design.mat =
      as.matrix(SNP.design.mat[, (k-1)*total.df.spline+1:total.df.spline])
    pred = as.matrix(cbind(nonSNP.design.mat, current.SNP.design.mat))
    resp = y.all
    for (i in 1:n.sub) {
      id.i = (id.all == i)
      W.0.5.i = W.0.5.all[(i-1)*n.rep + 1:n.rep, ]
      pred[id.i, ] = W.0.5.i %*% pred[id.i, ]
      resp[id.i] = W.0.5.i %*% resp[id.i]
    }
    resp.fitted = fitted(lm(resp ~ 0 + pred))
    u.vec[k] = mean((resp.fitted - resp)^2)
  }
  
  u.vec
}

## For this version, SNP.all does not have the same number of rows as others.
## Only 1 row in SNP.all is associated with each subject.

wanghuanMethod.2.parallel = function(y.all, time.all, gender.all, SNP.all.2, n.sub, n.rep, 
                            degree.spline = 3, knot.spline = 4, Mv = 15) {
  
  cluster = makeCluster(detectCores(), type = "SOCK")
  registerDoSNOW(cluster)
  
  # basic dimensions
  SNP.all.2 = as.matrix(SNP.all.2)
  p = ncol(SNP.all.2)
  id.all = rep(1:n.sub, each = n.rep)
  
  # generate b-spline functions for all the times
  total.df.spline = 1 + knot.spline + degree.spline
  time.spline = cbind(1, bs(time.all, df = total.df.spline - 1, degree = degree.spline, intercept = F))
  
  # generate the design matrix for non-SNP predictors
  nonSNP.repeated = cbind(1, gender.all)[, rep(1:2, each = total.df.spline)]
  time.spline.repeated = time.spline[, rep(1:total.df.spline, 2)]
  nonSNP.design.mat = nonSNP.repeated * time.spline.repeated
  
  # fit model (6) in the paper assuming i.i.d. error
  model6 = lm(y.all ~ 0 + nonSNP.design.mat)
  var.estimated.all = fitted(lm(model6$residuals^2 ~ 0 + time.spline)) # MAY BE NEGATIVE??
  
  # divide the estimated variances from model (6) so that the variances are all 1
  y.all.1 = y.all / sqrt(var.estimated.all)
  nonSNP.design.mat.1 = nonSNP.design.mat / sqrt(var.estimated.all)
  
  # fit model (6) in the paper assuming stationary M-dependent correlation structures,
  # the gee function requires that the numbers of observations in the groups are the same
  # check this!
  model6.2 = gee(y.all.1 ~ 0 + nonSNP.design.mat.1, id = id.all, family = gaussian("identity"), corstr ="non_stat_M_dep", Mv = Mv)
  cor.mat.estimated = model6.2$working.correlation
  cor.mat.estimated.inv = solve(cor.mat.estimated)
  
  # store the individual weight matrixes
  W.0.5.all = foreach(i=1:n.sub, .combine = rbind, .packages = "expm") %dopar% {
    id.i = (id.all == i)
    var.estimated.i = var.estimated.all[id.i]
    var.estimated.i.inv.5.mat = diag(1/sqrt(var.estimated.i))
    W.i = var.estimated.i.inv.5.mat %*% cor.mat.estimated.inv %*% var.estimated.i.inv.5.mat / n.rep
    sqrtm(W.i)
  }
  
  # fit model (5) in the paper for each kth SNP
  u.vec = foreach(k=1:p, .combine = c) %dopar% {
    
    SNP.i.design.mat = SNP.all.2[rep(1:n.sub, each = n.rep), rep(k, total.df.spline)] * time.spline
    pred = as.matrix(cbind(nonSNP.design.mat, SNP.i.design.mat))
    resp = y.all
    for (i in 1:n.sub) {
      id.i = (id.all == i)
      W.0.5.i = W.0.5.all[(i-1)*n.rep + 1:n.rep, ]
      pred[id.i, ] = W.0.5.i %*% pred[id.i, ]
      resp[id.i] = W.0.5.i %*% resp[id.i]
    }
    resp.fitted = fitted(lm(resp ~ 0 + pred))
    mean((resp.fitted - resp)^2)
  }
  
  stopCluster(cluster)
  u.vec
}





