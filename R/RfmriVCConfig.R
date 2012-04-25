#library(splines)
#library(Rniftilib)
#library(mgcv)
#library(mvtnorm)

#source("Z:/analyzes/Rtools/Rpack/RfmriVC-files/Fct_hrf-canonical_norm.r")
#source("Z:/analyzes/Rtools/Rpack/RfmriVC-files/Fct_designmatrix.r")
#source("Z:/analyzes/Rtools/Rpack/RfmriVC-files/Fct_bayes.r")
#source("Z:/analyzes/Rtools/Rpack/RfmriVC-files/Fct_ests_beta.R")
#source("Z:/analyzes/Rtools/Rpack/RfmriVC-files/Fct_ests_hrf.R")
#source("Z:/analyzes/Rtools/Rpack/RfmriVC-files/RfmriVCResults.r")

setOldClass(c("nifti"))

########################################################################################################
#define function validConfig for a validitity check of ConfigVC Objects
########################################################################################################

validConfig<-function(object){ 
	if(object@tr<=0)return("Error in ConfigVC: tr must be a positive number.")
	
	if(sum(dim(object@confounder))>2) if(dim(object@confounder)[1]!=dim(object@y)[4]) return("Error in ConfigVC: Row size of confounder matrix must match no. of fMRI images.")
	
	if(object@freqHPF<0)return("Error in ConfigVC: freqHPF must be a non-negative number.")
	
	
	if(object@nStimBF<1 || object@nStimBF>3){
		return("nStimBF must be 1, 2 or 3.")
	}
	
	if(length(object@stimTimes)!=length(object@VC)) return("Error in ConfigVC: For each stimulus presentation time point a value of the varying coefficient variable must be specified. Hence, length(stimTimes) must be equal to length(VC).")
	
	if(object@method!="classical" && object@method!="bayes") return(paste("Error in ConfigVC: Unknown method: ", object@method, ". Must either be 'classical' or 'bayes'."))
	
	return(TRUE)
}

########################################################################################################
# Class definition
########################################################################################################

setClass("ConfigVC", 
		representation(tr="numeric",
				stimTimes="vector",
				nStimBF="numeric",
				nVCknots="numeric",
				degreeVCspline="numeric",
				VC="vector",
				confounder="matrix",
				freqHPF="numeric",
				y="nifti",
				mask="nifti",	
				mask2totalIndex="matrix",
				total2maskIndex="nifti",
				grandMeanScaling="logical",
				U="matrix",
				X="matrix",
				penMat="matrix",
				method="character",
				methodParams="list"),
		validity=validConfig
)







#########################################################################################################
# Setter for common arguments and, by this, classical analysis (which does not need extra arguments)	#
#########################################################################################################

setGeneric("setConfigVC",function(object,value,...){standardGeneric("setConfigVC")})
setMethod(
		f="setConfigVC",
		signature="ConfigVC",
		function(object,tr=2,stimTimes,nStimBF=1,VC,nVCknots=10,degreeVCspline=2,confounder=NULL,freqHPF=128,y,mask=NULL,grandMeanScaling=TRUE){
			###################################################################################
			# Set arguments, variables an helper objects - if not specified set it to default #
			###################################################################################	
			object@tr<-tr
			object@stimTimes<-as.matrix(stimTimes)
			object@nStimBF<-nStimBF
			object@nVCknots<-nVCknots
			object@degreeVCspline<-degreeVCspline
			object@VC<-VC
			
			if(is.null(confounder)){
				object@confounder<-matrix(0)
			}else{
				object@confounder<-confounder
			}
			object@freqHPF<-freqHPF
			
			object@y<-y
			
			#Check if mask is specified, if not: generate a 'threshold' mask (default)
			if(is.null(mask)){
				mask<-nifti.image.new()
				mask$dim<-c(y$dim[1:3],1)   #same x,y,z- dimension as y
				nifti.image.alloc.data(mask)
				
				grandMeanDiv8<-mean(y[,,,,1,1,1])/8
				
				#Just select voxels with all time series values > grandMeanDiv8
				for(i in 1:dim(y)[1]){
					for(j in 1:dim(y)[2]){
						for(k in 1:dim(y)[3]){
							ifelse(sum(y[i,j,k,,1,1,1]<grandMeanDiv8)>0,mask[i,j,k,1,1,1,1]<-0,mask[i,j,k,1,1,1,1]<-1)
						}
					}
				}
			}	
			object@mask<-mask
			
			# Matrix with mapping from mask index to total index, 
			#	i.e. information about the i-th selected voxel is contained in row i, which contains the 3D index of the corresponding voxel
			m2tMat<-matrix(c(0,0,0),nrow=1,ncol=3)
			for(i in 1:dim(y)[1])for(j in 1:dim(y)[2])for(k in 1:dim(y)[3])if(mask[i,j,k,1,1,1,1]==1)m2tMat<-rbind(m2tMat,matrix(c(i,j,k),nrow=1,ncol=3))
			object@mask2totalIndex<-matrix(m2tMat[-1,],nrow=sum(mask[,,,1,1,1,1]==1),ncol=3)
			
			# Nifti image, which contains at i,j,k the index of the voxel within the set of selected voxels if selected by mask or zero otherwise
			total2maskIndex<-nifti.image.new()
			total2maskIndex$dim<-c(y$dim[1:3],1)   #same x,y,z- dimension as y
			nifti.image.alloc.data(total2maskIndex)
			voxCount<-1
			for(i in 1:dim(y)[1])for(j in 1:dim(y)[2])for(k in 1:dim(y)[3])if(mask[i,j,k,1,1,1,1]==1){
							if(mask[i,j,k,1,1,1,1]==1){
								total2maskIndex[i,j,k,1,1,1,1]<-voxCount
								voxCount<-voxCount+1
							}
						}
			object@total2maskIndex=total2maskIndex
							
			object@grandMeanScaling<-grandMeanScaling
			object@method<-"classical"
			object@methodParams<-list(NULL)
			
			####################
			# Prepare the data #
			####################	
			
			T<-dim(y)[4]
			t.grid<-1:T
			

			# Calculate Design matrix components 
			#------------------------------------
			
			#Calculate standardized HPF Baseline -> design matrix component U
			U<-as.matrix(w.design(t.grid,object@freqHPF,object@tr))
			#U<-cbind(U[,1],U[,-1]%*%diag(1/apply(U[,-1],2,sd)))
			#All DCT basis functions (except intercept) have same sd -> take sd of 2nd bf to standardize
			U<-U/sd(U[,2])
			#If available, attach standardized confounders to U 
			if(sum(dim(object@confounder))>2){
				C<-object@confounder
				Mean.C <- rep(1,T)%*%t(apply(C,2,mean))
				C<-(C-Mean.C)%*%diag(1/apply(C,2,sd))
				U<-cbind(U,C)
			}
			object@U<-U
			
			#Calculate design matrix part for varying stimulus coefficient
			X <- x.design(T=T, stimTimes=object@stimTimes, tr=object@tr, J.vig=object@VC, K=object@nVCknots, l=object@degreeVCspline,nStimBF=object@nStimBF)
			object@X<-X
			

			# Calculate penalty matrix 
			#---------------------------
			object@penMat<-k.pen(object@nVCknots+object@degreeVCspline-1,2)
			
			
			# Prepare y -> grand mean scaling 
			#----------------------------------
			#TODO: Do this before masking!?
			grandMean<-mean(object@y[,,,,1,1,1])
			object@y[,,,,1,1,1]<-100*object@y[,,,,1,1,1]/grandMean
				
			validObject(object)
			return(object)
		}
)


#############################################################################
# Setter for classical analysis												#			
# If needed, use this to reset ConfigVC-object to classical analysis		#
#############################################################################
setGeneric("setConfigVC.Classical",function(object,value,...){standardGeneric("setConfigVC.Classical")})
setMethod(
		f="setConfigVC.Classical",
		signature="ConfigVC",
		function(object){
			object@method<-"classical"
			object@methodParams<-list(NULL)
			validObject(object)
			return(object)
		}
)

###############################################################################
# Setter for bayesian analysis												  #	
# This setter changes the requested analysis type in object@method to "bayes" #
# If needed, non-default parameters for bayesian analysis can be assigned     # 
###############################################################################
setGeneric("setConfigVC.Bayes",function(object,value,...){standardGeneric("setConfigVC.Bayes")})
setMethod(
		f="setConfigVC.Bayes",
		signature="ConfigVC",
		function(object,chainSize=1100,burnin=100,step=5,sigma2start=1,sigma2IG=list(a=0.001,b=0.001),VCcontrol=NULL,diagSample=NULL){
			object@method<-"bayes"
			
			if(chainSize<=burnin){
				print("Error in setConfigVC.Bayes: chainSize must be larger than burnin. Parameters have been set to default sizes.")
				chainSize<-1100
				burnin<-100
				step<-5
			}
			
			remain<-(chainSize-burnin)/step
			if(remain<200){
				print(paste("Error in setConfigVC.Bayes: The value of 'chainsize' or 'burnin' or 'step' results in a chain size of", remain," for estimating each model parameter.  Parameters have been set to default sizes."))
				chainSize<-1100
				burnin<-100
				step<-5
			} 
			
			p<-object@nVCknots+object@degreeVCspline-1 #Number of gamma parameters per basis function
			element<-list(gammaStart=matrix(0,nrow=p,ncol=1),a=0.001,b=0.001) #default element of VCcontrol
			#If VCcontrol is null, assign default-element to each basis function specific list entry
			if(is.null(VCcontrol)|length(VCcontrol)!=object@nStimBF){
				if(!is.null(VCcontrol) & length(VCcontrol)!=object@nStimBF)print(paste("Warning: No. of VCcontrol elements",length(VCcontrol),"does not match nStimBF",object@nStimBF," -> reset to default."))
				VCcontrol<-list(element)
				k<-2
				while(k<=object@nStimBF){
					VCcontrol[[length(VCcontrol)+1]]<-element
					k<-k+1
				}
			}
			
			#If - for a basis function list element - an entry is falsely specified, replace it with default
			k<-1
			while(k<=object@nStimBF){
				if(prod(dim(VCcontrol[[k]][[1]]))!=p){
					print(paste("Warning: Length of starting VC coefficient vector for basis function",k,"has wrong dimension -> reset to default."))
					VCcontrol[[k]][[1]]=matrix(0,nrow=p,ncol=1)
				}
				if(VCcontrol[[k]][[2]]<=0){
					print(paste("Warning: IG shape parameter for basis function",k,"is not positive -> reset to default."))
					VCcontrol[[k]][[2]]=0.001
				}
				if(VCcontrol[[k]][[3]]<=0){
					print(paste("Warning: IG scale parameter for basis function",k,"is not positive -> reset to default."))
					VCcontrol[[k]][[3]]=0.001
				}
				k<-k+1
			}

			if(is.null(diagSample)){
				ifelse(sum(object@mask[,,,1,1,1,1])<20,nSamp<-sum(object@mask[,,,1,1,1,1]),nSamp<-20)
				diagSample<-sort(sample(1:sum(object@mask[,,,1,1,1,1]),nSamp))
			} 
						
			paramsList<-list(chainSize=chainSize,burnin=burnin,step=step,sigma2start=sigma2start,sigma2IG=sigma2IG,VCcontrol=VCcontrol,diagSample=diagSample)
			object@methodParams<-paramsList
			validObject(object)
			
			return(object)
		}
)

  






#########################################################################
# Function to conduct classical varying-coefficient fMRI analysis		#
#########################################################################
#setGeneric("doFmriVC.Classical",function(object,value,...){standardGeneric("doFmriVC.Classical")})
setGeneric("doFmriVC.Classical",function(object){standardGeneric("doFmriVC.Classical")})
setMethod(
		f="doFmriVC.Classical",
		signature="ConfigVC",
		function(object){
			if(!is.null(object@methodParams[[1]]))object<-setConfigVC.Classical(object)
			if(object@method!="classical")object<-setConfigVC.Classical(object)
			print("Classical varying-coefficient fMRI analysis")
			dimU<-dim(object@U)
			dimX<-dim(object@X)
			#complete design matrix
			X.total<-cbind(object@U,object@X)
			
			#Penalty parameter
			if(object@nStimBF==3){
				SpenMat<-list(object@penMat,object@penMat,object@penMat)
				off<-c(dimU[2]+1,dimU[2]+dimX[2]/3+1,dimU[2]+dimX[2]*2/3+1)
				sp<-c(1,1,1)
			}else{
				if(object@nStimBF==2){
					SpenMat<-list(object@penMat,object@penMat)
					off<-c(dimU[2]+1,dimU[2]+dimX[2]/2+1)
					sp<-c(1,1)
				}else{
					SpenMat<-list(object@penMat)
					off<-dimU[2]+1
					sp<-1
				}
			}
			
			#Initialize nifit objects for estimates 
			thetaHat<-nifti.image.new()
			thetaHat$dim<-c(object@y$dim[1:3],dim(X.total)[2]) 
			nifti.image.alloc.data(thetaHat)
			
			sigma2Hat<-nifti.image.new()
			sigma2Hat$dim<-c(object@y$dim[1:3],1)   #same x,y,z- dimension as y
			nifti.image.alloc.data(sigma2Hat)
			
			lambdaHat<-nifti.image.new()
			lambdaHat$dim<-c(object@y$dim[1:3],object@nStimBF)   #same x,y,z- dimension as y
			nifti.image.alloc.data(lambdaHat)
			
			voxCount<-1
			dims<-dim(object@y)
			print(paste(sum(object@mask[,,,1,1,1,1] == 1),"voxels selected by mask for analysis of a total of",prod(dims[1:3]),"voxels."))
			for(i in 1:dims[1]){
				
				if(voxCount<=sum(object@mask[,,,1,1,1,1] == 1))print(paste("Processing voxel",voxCount,":",Sys.time()))
				
				for(j in 1:dims[2]){
					for(k in 1:dims[3]){
						
						# Do Voxelwise estimation if voxel is selected by mask
						if(object@mask[i,j,k,1,1,1,1] == 1){
							# Voxelwise estimation
					
							vcmHat <- magic(y=object@y[i,j,k,,1,1,1],X=X.total,sp=sp,S=SpenMat,off=off)
							#vcmHat <- mgcv(y=object@y[i,j,k,,1,1,1],X=X.total,sp=sp,S=SpenMat,off=off)
					
				
							# save beta
							thetaHat[i,j,k,,1,1,1] <- vcmHat$b
							
							# save sigma^2
							sigma2Hat[i,j,k,1,1,1,1] <- vcmHat$scale
							
							# save lambda 
							lambdaHat[i,j,k,,1,1,1] <- vcmHat$sp
							
							voxCount<-voxCount+1
						}
					}
				}
			}
			print(paste("Classical algorithm finished:",Sys.time()))
			resultsObj<-new("ResultsVC")
			resultsObj<-setResultsVC(object=resultsObj,configVCobj=object,method="classical",methodResults=list(thetaHat=thetaHat,sigma2Hat=sigma2Hat,lambdaHat=lambdaHat))
			return(resultsObj)
		}
)

#########################################################################
# Function to conduct bayesian varying-coefficient fMRI analysis		#
#########################################################################
#setGeneric("doFmriVC.Bayes",function(object,value,...){standardGeneric("doFmriVC.Bayes")})
setGeneric("doFmriVC.Bayes",function(object){standardGeneric("doFmriVC.Bayes")})
setMethod(
		f="doFmriVC.Bayes",
		signature="ConfigVC",
		function(object){
			if(is.null(object@methodParams[[1]]))object<-setConfigVC.Bayes(object)
			if(object@method!="bayes")object<-setConfigVC.Bayes(object)
			
			print("Bayesian varying-coefficient fMRI analysis")
			nVoxSel<-sum(object@mask[,,,1,1,1,1] == 1) #No. of voxels selected by mask for analysis
			p<-length(object@methodParams$VCcontrol[[1]][[1]])[1] #Length of a singel, voxelwise gamma vector
			remain<-(object@methodParams$chainSize-object@methodParams$burnin)/object@methodParams$step #No. of iterations used for estimation
			nStimBF<-object@nStimBF 
			
			#Array for storing gamma trajectories (for slected voxels only)
			gammaSamples<-array(dim=c(nVoxSel,p,remain,nStimBF))
			
			#Initialize nifti objects for thetaHat estimates for c(zeta,gammaVec) 
			thetaHat<-nifti.image.new()
			thetaHat$dim<-c(object@y$dim[1:3],dim(object@U)[2]+nStimBF*p) 
			nifti.image.alloc.data(thetaHat)
			
			#TODO: Wieder entfernen?
			diagnostics<-list(NULL)	#list for saving parameter trajectories
			d<-1					#Itertes through voxels selected for diagnostics
			
			voxCount<-1				#Iterates through voxels selected for analysis 
			dims<-dim(object@y)
			print(paste(nVoxSel,"voxels selected by mask for analysis of a total of",prod(dims[1:3]),"voxels."))
			for(i in 1:dims[1]){
				if(voxCount<=sum(object@mask[,,,1,1,1,1] == 1))print(paste("Processing voxel",voxCount,":",Sys.time()))
				for(j in 1:dims[2]){
					for(k in 1:dims[3]){
						
						# Do Voxelwise estimation if voxel is selected by mask
						if(object@mask[i,j,k,1,1,1,1] == 1){				
							gibbs <-VCM.gibbs(y=object@y[i,j,k,,1,1,1],object@U,object@X,object@penMat,object@methodParams)
							
							
							###########################
							# Store results
							###########################				
							
							#Store gamma trajectory in array for by mask selected voxels 
							gammaSamples[voxCount,,,]<-gibbs$gamma
							gammaHat<-apply(gibbs$gamma,c(1,3),mean)
							zetaHat<-apply(gibbs$zeta,1,mean)
							thetaHat[i,j,k,,1,1,1]<-c(zetaHat,as.vector(gammaHat))
							
							
							#Store trajectories if voxel is selected for diagnostics
							if(sum(voxCount==object@methodParams$diagSample)){
								diagnostics[[d]]<-gibbs
								d<-d+1
							}
							
							voxCount<-voxCount+1
						}
					}
				}
			}
			print(paste("MCMC algorithm finished:",Sys.time()))
			resultsObj<-new("ResultsVC")
			resultsObj<-setResultsVC(object=resultsObj,configVCobj=object,method="bayes",methodResults=list(thetaHat=thetaHat,gammaSamples=gammaSamples,diagnosticsSel=diagnostics))
			return(resultsObj)
		}
)

#############################################################
# Function to conduct varying-coefficient fMRI analysis		#
# Automatic selection of classical or bayesian analysis		#
#############################################################
#setGeneric("doFmriVC",function(object,value,...){standardGeneric("doFmriVC")})
setGeneric("doFmriVC",function(object){standardGeneric("doFmriVC")})
setMethod(
		f="doFmriVC",
		signature="ConfigVC",
		function(object){
			
			if(object@method=="bayes"){
				resultsObj<-doFmriVC.Bayes(object)
			}else{
				resultsObj<-doFmriVC.Classical(object)
			}
			
			return(resultsObj)
		}
)


########################################################################################################
#Show
setMethod(
		f="show",
		signature="ConfigVC",
		function(object){
			cat("This is an object of Class 'ConfigVC'. Package 'RfmriVC'.\n")
			cat("It contains the following configuration options for a call to the RfmriVC-algorithm.\n\n")
			
			cat(paste("method = Requested algorithm:",object@method,"\n\n"))
			
			cat(paste("tr = Repetition time of fMRI scans:",object@tr,"\n"))
			cat(paste("stimTimes = Vector of stimulus presentation times: "))
			cat(paste(head(object@stimTimes)))
			cat(" ...\n")
			cat(paste("nStimBF = Number of canonical HRF Basis functions:",object@nStimBF,"\n"))
			cat(paste("For each effect of the ",object@nStimBF," HRF basis functions a VC is requested.\n"))
			cat(paste("VC = Variable the varying coefficient depends on: "))		
			cat(paste(head(object@VC)))
			cat(" ...\n")
			cat(paste("nVCknots = Number of knots for spline basis of VC:",object@nVCknots,"\n"))
			cat(paste("degreeVCspline = Degree of the VC spline:",object@degreeVCspline,"\n"))
			cat(paste("Design matrix part for varying stimulus coefficients are contained in slot X.\n"))
			cat(paste("confounder = confounder matrix with dimension:",dim(object@confounder)[1]," , ",dim(object@confounder)[2],"\n"))
			cat(paste("freqHPF = Frequency of DCT high pass filter:",object@freqHPF,"\n"))
			cat(paste("Nuisance regressors (HPF and confounders) are contained in slot U.\n"))
			
			cat(paste("\ny = 4D nifti file with fMRI data: \n"))
			cat(paste(show(object@y)))
			cat(paste("\nmask = 3D nifti file with analysing mask:\n"))
			cat(paste(show(object@mask)))
			cat(paste("\nGrand mean scaling:",object@grandMeanScaling,"\n"))
			
			cat("Further slots: penMat (penalty matrix for one VC spline), mask2totalIndey (mapping from mask to total index), total2maskIndex (mapping from total to mask index)\n")
			
			
			if(object@method=="bayes"){
				cat(paste("Configuration parameters for bayesian analysis: \n"))
				cat(paste("MCMC chainsize: ",object@methodParams$chainSize,"\n"))
				cat(paste("MCMC burnin: ",object@methodParams$burnin,"\n"))
				cat(paste("MCMC thinning: ",object@methodParams$step,"\n"))
				cat(paste("MCMC sigma2 start: ",object@methodParams$sigma2start,"\n"))
				cat(paste("MCMC sigma2 IG parameters: a = ",object@methodParams$sigma2IG[[1]], ", b = ",object@methodParams$sigma2IG[[2]],"\n"))
				cat(paste("MCMC VC control BF",c(1,2,3),": ",object@methodParams$VCcontrol,"\n"))
			}
		}
)


