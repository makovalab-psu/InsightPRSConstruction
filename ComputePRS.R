#################################################################################
# Functions to compute PRS from FLAME results
#################################################################################

#--------------------------------------------------------------------------------
# load necessary libraries
#--------------------------------------------------------------------------------
library(fda)
library(CompQuadForm)
library(fdapace)
source("FS_penreg_v2.R")

#--------------------------------------------------------------------------------
# function to compute pairwise integrals for finding weights
#--------------------------------------------------------------------------------
compute_pairwise_integrals <- function(matrix, T_grid)
{
    matrix_derivatives <- apply(matrix, 2, function(y)
    {
        t(as.vector(apply(matrix,2,
                          function(x)
                          {
                              sum(x*y*(T_grid[2]-T_grid[1]))
                          }
        )
        ))
    }
    )

    return(matrix_derivatives)
}
l2 <- function(v){
  return(sqrt(sum(v^2)))
}


#--------------------------------------------------------------------------------
# function to construct the PRS and evaulate on curves 
# input:
# Y.f - functional object representing curves
# X - nxp design matrix (columns would be SNPs)
# mod.flame - object with final FLAME model components
# output:
# weights - weights for each selected variable (SNPs)
# PRS - score assigned to each individual (nx1 vector)
# pvalues - p-values from regressing PRS on Y.f. Evaluated 3 ways
# R2overtime - R^2 "curve" across each time point
# R2all - overall R^2
#--------------------------------------------------------------------------------
ConstructPRS <- function(Y.f, X, mod.flame){
  idx_pred = mod.flame$predictors
  beta = mod.flame$beta
  Xflame = X[,idx_pred]

  # get original time domain from functional object
  T_domain = sort(c(Y.f$basis$rangeval, Y.f$basis$params))

  # Evaluate at time points to get functions (rather then projected values)
  Bt = eval.fd(T_domain,beta) # returns mxp matrix (zeros where variable not selected)
  Bt = Bt[,idx_pred] # can cut out other variables

  # get sample covariance
  # change estimation method if n << length(idx_pred)
  A = cov(Xflame)

  # Compute integral of pairwise product of functions
  # Note: function taken from flame package
  B = compute_pairwise_integrals(Bt,T_domain)
  C = A%*%B%*%A

  # Maximize objective t(w)CW subject to t(w)*w = 1
  # Solution is eigenvector corresponding to max eigen value since we have a quadratic form
  # and C,I are symmetric
  # Note: accurate up to a sign
  w_cov = -eigen(C)$vectors[,1]

  # Weight and sum the values using w_cov
  W_cov = matrix(,nrow=dim(Xflame)[1],ncol=dim(Xflame)[2])
  for(i in 1:dim(Xflame)[2]){
          W_cov[,i] = w_cov[i]*Xflame[,i]
  }
  PRS_cov = rowSums(W_cov)

  # add on column for intercept and create X matrix
  # Note: intercept not really necessary if Y.f has been centered already
  Xmat_cov = as.matrix(cbind(1,PRS_cov))
  # function-on-scalar regression
  myfit_cov = FS_penreg(Y.f,Xmat_cov,q=20)

  # Calculate p-values and get coefficients
  pv_cov = myfit_cov$pvals[-1,c(3,5,7)]
  Bhat_f_cov = myfit_cov$Bhat[2]
  se_fun_cov = myfit_cov$se_fun[2]
  low_cov = Bhat_f_cov - 2*se_fun_cov
  high_cov = Bhat_f_cov + 2*se_fun_cov
  myfuns_cov = Bhat_f_cov
  myfuns_cov$coefs = cbind(Bhat_f_cov$coefs, low_cov$coefs, high_cov$coefs)

  # R^2 measurement across time domain
  Yt = eval.fd(T_domain,Y.f)
  Bt_cov = eval.fd(T_domain,myfit_cov$Bhat)
  Bt_cov_se = eval.fd(T_domain,myfit_cov$se_fun)
  Ypred_cov = Xmat_cov%*%t(Bt_cov)
  R2 = vector(,length = length(T_domain))
  for(i in 1:length(T_domain)){
    R2[i] = 1-(sum((Ypred_cov[,i]-Yt[i,])^2))/sum((Yt[i,]-mean(Yt[i,]))^2)
  }
  
  # get one R^2 measurement
  Yt = t(Yt)
  D = vector(,length=dim(Yt)[1])
  D2 = vector(,length=dim(Yt)[1])
  CMY = colMeans(Yt)
  for(i in 1:length(D)){
    D[i] = l2(Yt[i,]-Ypred_cov[i,])
    D2[i] = l2(Yt[i,]-CMY)
  }
  R2.1 = 1-(sum(D^2))/(sum(D2^2))

  return(list(weights = w_cov,
              PRS = PRS_cov,
              pvalues = pv_cov,
              R2overtime = R2,
              R2all = R2.1))
}
