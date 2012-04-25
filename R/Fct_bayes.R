################################################################################
# Functions to conduct Bayesian Gibbs Sampling 
# (i.e. full conditional distribution samplers and whole Gibbs sampler)
# Created on: 14.11.11
# Last change: 24.04.12
# Author: Ludwig Bothmann
# Changed by: Stefanie Kalus
################################################################################

# random samples from full conditional distribution of fixed effects zeta
fc.zeta <- function(sigma2,gammaVec,y,X,U,C.inv,c){
	
	Sigma.zeta <- sigma2 * solve(t(U)%*%U + sigma2*C.inv)
	
	mu.zeta <- (1/sigma2) * Sigma.zeta %*% (t(U) %*% (y - X %*% gammaVec) + sigma2 * (C.inv %*% c))
	
	zeta.neu <- matrix(rmvnorm(1, mu.zeta, Sigma.zeta,method="svd"),nrow=dim(U)[2],ncol=1)
	
	return(zeta.neu)
}

# random samples from full conditional distribution of VC gamma.k (basis function component k)
fc.gamma.k <- function(y,gamma.k,X.k,tau2.k,sigma2,K.pen,eta){
	
	Sigma.gamma.k <- solve((1/sigma2) * t(X.k)%*%X.k + (1/tau2.k) * K.pen)
	
	eta.k <- eta - X.k%*%gamma.k
	
	mu.gamma.k <- (1/sigma2) * Sigma.gamma.k %*% (t(X.k) %*% (y - eta.k))
	
	gamma.k.neu <- rmvnorm(1, mu.gamma.k, Sigma.gamma.k)
	
	return(gamma.k.neu)
}

# random sample from full conditional distribution of sigma2
fc.sigma2 <- function(zeta,gammaVec,y,X,U,a_sigma,b_sigma){
	
	a.sigma2 <- a_sigma + (dim(X)[1]/2)

	eta1<-U %*% zeta
	eta2<-X %*% gammaVec
	eta <- eta1 + eta2  #Each summand has to be calculated separately, otherwise it does not work
	
	b.sigma2 <- b_sigma + (1/2) * sum((y-eta)^2)
	
	sigma2.neu <- 1/rgamma(1, shape = a.sigma2, rate = b.sigma2)
	
	return(sigma2.neu)
}

# random sample from full conditional distribution of tau2.k (basis function component k)
fc.tau2.k <- function(gamma.k,a.k,b.k,K.pen){
	
	a.tau2 <- a.k + qr(K.pen)$rank/2
	
	b.tau2 <- b.k + (1/2) * t(gamma.k) %*% K.pen %*% gamma.k
	
	tau2.neu <- 1/rgamma(1, shape = a.tau2, rate = b.tau2)
	
	return(tau2.neu)
}





################################################################################
# Gibbs-Sampler
#	- y		Vector with fMRI signal time series of one voxel
#   - U		design matrix of fixed effects (baseline+confounders)
#	- X		total design matrix of varying coefficients
#	- K.pen				Smoother matrix for a gamma.k update
#	- method.Params		Parameters specific to full bayesian analysis, like hyperparameters,
#						starting values, MCMC algorithm parameters 
#
################################################################################


VCM.gibbs <- function(y,U,X,K.pen,methodParams){
	nStimBF<-length(methodParams$VCcontrol) #nStimBF
	p<-length(methodParams$VCcontrol[[1]][[1]])[1] #Length of a gamma vector
	remain<-(methodParams$chainSize-methodParams$burnin)/methodParams$step #No. of iterations used for estimation
	
	#######################################################################
	# 1. Initialize matrices/objects for (non-discarded) samples
	#######################################################################
	
	# zeta
	zeta <- matrix(nrow=dim(U)[2],ncol=remain)
	
	# array for gamma=(gamma.k) with 1<=k<=nStimBF
	gamma<-array(dim=c(p,remain,nStimBF))
	
	# sigma^2
	sigma2 <- vector(length=remain)
	
	# matrix for tau2=(tau2.k) with 1<=k<=nStimBF
	tau2<-matrix(nrow=nStimBF,ncol=remain)

	
	#######################################################################
	# 2. Initialize objects which store actual state of parameters
	#######################################################################
	
	# matrix for gamma.l=(gamma.l.k) state at iteration l with 1<=k<=nStimBF
	gamma.l<-matrix(nrow=p,ncol=nStimBF)
	for(k in 1:nStimBF){
		gamma.l[,k]<-methodParams$VCcontrol[[k]][[1]]
	}
	
	# vector for tau2.l=(tau2.l.k) state at iteration l with 1<=k<=nStimBF
	tau2.l<-rep(NA,nStimBF)

	
	# vector for zeta.l state at iteration l
	zeta.l<-matrix(NA,nrow=dim(U)[2],ncol=1)
	
	# sigma2 state at iteration l
	sigma2.l <-methodParams$sigma2start
	


	
	#######################################################################
	# 3. Gibbs sampling 
	#######################################################################
	
	# Set prior matrices 
	# C.inv: zero matrix
	C.inv <- matrix(nrow=dim(U)[2],ncol=dim(U)[2],data=0)
	# c: vector with ones

	c <- rep(1,dim(U)[2])
	selIter<-1
	for(l in 1:methodParams$chainSize){		
		# tau^2_k for k = 1,2,3
		for(k in 1:nStimBF)tau2.l[k]<-fc.tau2.k(gamma.l[,k],methodParams$VCcontrol[[k]][[2]],methodParams$VCcontrol[[k]][[3]],K.pen)
		
		
		# zeta
		zeta.l<-fc.zeta(sigma2.l,as.vector(gamma.l),y,X,U,C.inv,c)
		
		# sigma^2
		sigma2.l<- fc.sigma2(zeta.l, as.vector(gamma.l), y,X,U,methodParams$sigma2IG$a,methodParams$sigma2IG$a)
		
		# gamma_k for k = 1,2,3
		for(k in 1:nStimBF){
			X.k<-X[,(1:p)+(k-1)*p]
			gammaVec<-as.vector(gamma.l)
			eta1<-U%*%zeta.l
			eta2<-X%*%gammaVec
			eta<-eta1+eta2 #Each summand has to be calculated separately, otherwise it does not work
			gamma.l[,k]<-fc.gamma.k(y,gamma.l[,k],X.k,tau2.l[k],sigma2.l,K.pen,eta)	
		}
		
		#Store only those states which are selected according to burnin and step parameter
		stepSel<-(l-methodParams$burnin-1)%%methodParams$step;
		if((l>methodParams$burnin) & (!stepSel)){
			sigma2[selIter]<-sigma2.l
			zeta[,selIter]<-zeta.l
			tau2[,selIter]<-tau2.l
			#print(gamma.l)
			#print(gamma[,selIter,])
			gamma[,selIter,]<-gamma.l
			#print(gamma[,selIter,])
			selIter<-selIter+1
		}
	}
	
	#######################################################################
	# 4. Return
	#######################################################################

	samples <- list(zeta = zeta, gamma = gamma, sigma2 = sigma2, tau2= tau2)	
	return(samples)
}



