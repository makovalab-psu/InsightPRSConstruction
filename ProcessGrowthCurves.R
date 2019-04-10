#################################################################################
# Functions to process growth curves
#################################################################################

#--------------------------------------------------------------------------------
# load necessary libraries
#--------------------------------------------------------------------------------
library(fda)
library(fdapace)


#--------------------------------------------------------------------------------
# Function to convert curves to functional object for FLAME and PRS construction
# input:
# Y - nxt matrix, each row corresponds to an individual with values over t points
# tm - nxt marix, each row corresponds to the time point where Y value is observed
# output:
# Y.f - functional object
#--------------------------------------------------------------------------------
ConvertFunction <- function(Y,tm){
	N = dim(Y)[1];M=dim(Y)[2]
	Lt = list(); Ly=list()
	for(i in 1:N){Lt[[i]] = tm[i,];Ly[[i]] = Y[i,]}
	FPCA_fit=FPCA(Ly,Lt,optns=list(userBwCov=5,userBwMu=5))
	Y_imp=fitted(FPCA_fit)
	Y_grid=FPCA_fit$workGrid
	# Center response
	Y_imp_c = scale(Y_imp,scale=FALSE)
	Y.f = Data2fd(Y_grid,t(Y_imp_c))
	return(Y.f)
}

#--------------------------------------------------------------------------------
# Function to convert growth values to the "shift" regime switching from W/H to
# BMI smoothly from ages 24 months to 36 months
# input:
# weight - weight at age where response is adjusted 
# height - height at age where response is adjusted 
# age - age (in months) where response is adjusted
# output:
# new - one value for new response at age
#--------------------------------------------------------------------------------
changeR <- function(weight,height,age){
  f = exp(-(1/((1/12)*age-2)))/(exp(-(1/((1/12)*age-2)))+exp(-(1/(1-((1/12)*age-2)))))
  new = weight/(height^(1+f))
  return(new)
}