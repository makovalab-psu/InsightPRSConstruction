#################################################################################
# Function to apply FLAME for variable selection
#################################################################################

#--------------------------------------------------------------------------------
# load necessary libraries
#--------------------------------------------------------------------------------
library(flm) # may need to specify local library installed
library(fda)
library(fdapace)

#--------------------------------------------------------------------------------
# Function to apply FLAME to data
# input:
# Y.f - functional object for observations
# X - nxp design matrix, each column would correspond to a SNP
# k - kill switch parameter in FLAME to limit maximum # of variables selected
# r - ratio to compute minimum value of lambda, default is 0.01 (used 0.01 in results)
# output:
# mod.flame - object containing final FLAME model estimates
# varnames.select - list of variables selected (names of SNPs)
# comptime - computing time for FLAME
#--------------------------------------------------------------------------------
applyFLAME <- function(Y.f,X,k,r){
	X = scale(X)
	ptm = proc.time()
	mod.flame = FLAME(Y.f, X=X, type_kernel='sobolev',number_non_zeros=k, ratio_lambda = r,verbose=TRUE)
	comptime = proc.time()-ptm
	idx_pred = mod.flame$predictors
	varnames.select = colnames(X)[idx_pred]
	return(list(mod.flame = mod.flame,
				varnames.select = varnames.select,
				comptime = comptime
				)
	)
}