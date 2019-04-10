######################################################################
# FS_penreg function is for function-on-scalar regression
# Used to regress PRS on curves
######################################################################


# This code carries out Function-on-Scalar regression.
# It does this using FPCA as a basis for the functions
# as well as penalizing the second derivative.

## Inputs ##
# Y.f: Functional Object
# Xmat: matrix of scalar predictors (N x R matrix)
# penlamb_all: vector of penalty parameters to consider.  
# The penalty parameter for each beta is selected by choosing
# the best parameter among penlamb_all using gcv. Thus each
# beta can have a different penalty parameter.  
# q: number of FPCs to use


## Outputs ##
# Bhat: coefficient matrix for functional beta that takes same basis as Y's
#		(R x K matrix)
# se_fun: the point-wise standard error functions for each beta
# pvals: test statistics size and p-values based on L2-test, PC-test, 
#		and Choi-test
# lambda: the lambdas that were chosen


FS_penreg <- function(Y.f, Xmat, penlamb_all=c(0,10^(0:10)),q=10) {
	#penlamb_all<-c(0,10^(0:10))
	#q=10
	Y.pc<-pca.fd(Y.f,nharm=q)
	v<-Y.pc$harmonics
	Ymat<-inprod(Y.f,v)
	v2 = deriv.fd(v,2)
	penmat = inprod(v2,v2)
  
  
  
  	N = dim(Ymat)[1]
  	K = dim(Ymat)[2]
  	R = dim(Xmat)[2]

  	if (dim(Xmat)[1] != N) {
    	stop("Check the dimension of Ymat and Xmat")
  	}
  	if (!is.matrix(Xmat)) {
    	stop("Xmat needs to be a matrix")
  	}
  	
  	penlamb<-rep(0,times=R)
  	Lamb <- diag(penlamb)
  	# Penalized Estimation
  	V = c(t(Ymat)%*%Xmat)
  	A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
  	Ainv = chol2inv(chol(A))
  	Bhat_v = Ainv%*%V
  	Bhat_m_T = matrix(Bhat_v,nrow=K,ncol=R)
  	Bhat = t(Bhat_m_T)	### Penalized LSE
  	Eps_mat = Ymat - Xmat%*%Bhat
  	df<-N-sum(diag(Ainv%*%((t(Xmat)%*%Xmat)%x%diag(1,nrow=q))))/K
  	GCV<-N*sum((Eps_mat^2))/df^2
  	
  	max_iter<-10;cnt<-0
	lamb_ind<-rep(1,times=R)
	go<-TRUE
  	while(go==TRUE & cnt<max_iter){
  		go<-FALSE
  		for(j in 1:R){
  			if(lamb_ind[j]<length(penlamb_all)){
	  			penlamb_new<-penlamb
	  			penlamb_new[j]<-penlamb_all[lamb_ind[j]+1]
	  			Lamb <- diag(penlamb)
  				# Penalized Estimation
  				V = c(t(Ymat)%*%Xmat)
  				A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
  				Ainv = chol2inv(chol(A))
  				Bhat_v = Ainv%*%V
  				Bhat_m_T = matrix(Bhat_v,nrow=K,ncol=R)
  				Bhat = t(Bhat_m_T)	### Penalized LSE
  				Eps_mat = Ymat - Xmat%*%Bhat
  				df<-N-sum(diag(Ainv%*%((t(Xmat)%*%Xmat)%x%diag(1,nrow=q))))/K
  				GCV_new<-N*sum((Eps_mat^2))/df^2
  				if(GCV_new<=GCV){
  					GCV<-GCV_new
  					penlamb<-penlamb_new
  					lamb_ind[j]<-lamb_ind[j]+1
  					go<-TRUE
  				}
			}
  		}
  		cnt<-cnt+1
  		
  	}
  	
  	
  	
  	
  	
  	
  	Lamb <- diag(penlamb)
  	# Penalized Estimation
  	V = c(t(Ymat)%*%Xmat)
  	A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
  	Ainv = chol2inv(chol(A))
  	Bhat_v = Ainv%*%V
  	Bhat_m_T = matrix(Bhat_v,nrow=K,ncol=R)
  	Bhat = t(Bhat_m_T)	### Penalized LSE
  	Eps_mat = Ymat - Xmat%*%Bhat
  	Sig_mat = (1/(N - R)) * t(Eps_mat)%*%Eps_mat
  	Big_Cov_mat = Ainv%*%((t(Xmat)%*%Xmat)%x%Sig_mat)%*%Ainv


  	# Testing
  	EV_PC<-.85
  	EV_Choi<-.99
  	EV_L2<-.999

  	pvalmat <- data.frame(betanum=1:nrow(Bhat), Ti_L2 = NA, pval_L2 = NA, 
						Ti_PC = NA, pval_PC = NA,
						Ti_Choi = NA, pval_Choi = NA)
  	require('CompQuadForm')
  	se_vec<-matrix(ncol=nrow(Bhat),nrow=50)
  	for (i in 1:nrow(Bhat)) {

    	covBi = Big_Cov_mat[(K*(i-1)+1):(K*i),(K*(i-1)+1):(K*i)]
    	Eig = eigen(covBi, symmetric=T)
    	eigvals_i <- Eig$values
    	eigvecs_i <- Eig$vectors

    	Scores<-Bhat%*%eigvecs_i

    	k_PC<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_PC)
    	k_Choi<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_Choi)
    	k_L2<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_L2)
    	
  
  
    	Ti_L2 = sum( Scores[i,1:k_L2]^2 )
    	Ti_PC = sum( Scores[i,1:k_PC]^2/eigvals_i[1:k_PC] )
    	Ti_Choi = sum( Scores[i,1:k_Choi]^2/sqrt(eigvals_i[1:k_Choi]) )
    	pvalmat$Ti_L2[i] = Ti_L2
    	pvalmat$Ti_PC[i] = Ti_PC
    	pvalmat$Ti_Choi[i] = Ti_Choi
    	pvalmat$pval_L2[i] = imhof(Ti_L2/eigvals_i[1], eigvals_i[1:k_L2]/eigvals_i[1])$Qq
    	pvalmat$pval_PC[i] = pchisq(Ti_PC, df = k_PC, lower.tail=FALSE)
    	pvalmat$pval_Choi[i] = imhof(Ti_Choi/sqrt(eigvals_i[1]), sqrt(eigvals_i[1:k_Choi])/sqrt(eigvals_i[1]),epsabs=10^(-16))$Qq
      
    	Cv_c<- v$coefs%*%covBi%*%t(v$coefs)
    	Cv_fun<-bifd(Cv_c,v$basis,v$basis)
    	pts<-seq(v$basis$rangeval[1],v$basis$rangeval[2],length=50)
    	Cv_mat = eval.bifd(pts,pts,Cv_fun)
    	se_vec[,i]<-sqrt(diag(Cv_mat))
    	
  	}
  	Bhat_f<-fd(coef=v$coefs%*%t(Bhat),Y.f$basis)
  	se_fun<-Data2fd(pts,se_vec,Y.f$basis)
  	return(list(Bhat = Bhat_f, se_fun=se_fun, pvals = pvalmat,lambda=penlamb))
}

