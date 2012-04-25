
########################################################################################################
#define function validResults for a validitity check of ResultsVC Objects
########################################################################################################

validResults<-function(object){ 
	if(object@alpha<=0 | object@alpha>=1)return("Error in ResultsVC: signifiance level alpha must be a value in (0,1)")
	
	if(length(object@stimTimes)!=length(object@VC)) return("For each stimulus presentation time point a value of the varying coefficient variable must be specified. Hence, length(stimTimes) must be equal to length(VC).")
	
	if(object@method!="classical" && object@method!="bayes") return(paste("Unknown method: ", object@method, ". Must either be 'classical' or 'bayes'."))
	
	if(min(object@VC.plot)<min(object@VC)) return(paste("VC.plot values are out of bounds: ", min(object@VC.plot)," < ", min(object@VC) ))
	if(max(object@VC.plot)>max(object@VC)) return(paste("VC.plot values are out of bounds: ", max(object@VC.plot)," > ", max(object@VC) ))
	
	return(TRUE)
}

########################################################################################################
# Class definition
########################################################################################################

setClass("ResultsVC", 
		representation(VC="vector",
				VC.plot="vector",
				nVCknots="numeric",
				degreeVCspline="numeric",
				stimTimes="vector",
				nStimBF="numeric",
				U="matrix",
				X="matrix",
				penMat="matrix",
				total2maskIndex="nifti",
				t.plot="vector",
				alpha="numeric",
				method="character",
				methodResults="list"),
		validity=validResults
)


#########################################################################################################
# Setter for parameters needed for post-hoc analysis - used as return of doFmriVC-function				#
#########################################################################################################

setGeneric("setResultsVC",function(object,value,...){standardGeneric("setResultsVC")})
setMethod(
		f="setResultsVC",
		signature="ResultsVC",
		function(object,configVCobj,method,methodResults){
			object@VC=configVCobj@VC
			object@VC.plot=seq(min(configVCobj@VC),max(configVCobj@VC),length.out=1000)
			object@nVCknots=configVCobj@nVCknots
			object@degreeVCspline=configVCobj@degreeVCspline
			object@stimTimes=configVCobj@stimTimes
			object@nStimBF=configVCobj@nStimBF
			object@U=configVCobj@U
			object@X=configVCobj@X
			object@penMat=configVCobj@penMat
			object@total2maskIndex=configVCobj@total2maskIndex
			object@t.plot=seq(0,30,length.out=1000)
			object@alpha=0.05
			object@method=method
			object@methodResults=methodResults
			validObject(object)
			return(object)
		}
)

#########################################################################################
# Parameters whose settings can be customized by user can be reset with this function.	#
# The remainder parameters are set directly in the analysing routine  					#
#########################################################################################

setGeneric("resetResultsVC",function(object,value,...){standardGeneric("resetResultsVC")})
setMethod(
		f="resetResultsVC",
		signature="ResultsVC",
		function(object,VC.plot=seq(min(configVCobj@VC),max(configVCobj@VC),length.out=1000),t.plot=seq(0,30,length.out=1000),alpha=0.05){
			object@VC.plot=sort(VC.plot)
			object@t.plot=sort(t.plot) 
			object@alpha=alpha
			validObject(object)
			return(object)
		}
)



#############################################################################################################################
# Post-hoc calculations																										#
#############################################################################################################################

#############################################################################
# Function to calculate point and interval estimates of beta_i(VC) for a
# specified grid VC.plot for a single voxel i
# 
# Arguments: 
#	- object		object of class ResultsVC
#	- voxelIndex3D	3D vector with 3D index of requested voxel 
# 	- VC.plot		(optional) grid of VC variable values for which betaHat should be predicted					
#	- which.beta	corresponding HRF basis function number (default: canonical HRF <=> which.beta=1)
#############################################################################

setGeneric("betaHat",function(object,value,...){standardGeneric("betaHat")})
setMethod(
		f="betaHat",
		signature="ResultsVC",
		function(object,voxelIndex3D,VC.plot=object@VC.plot,which.beta=1,alpha=object@alpha){
			
			p<-object@nVCknots+object@degreeVCspline-1 #Number of gamma parameters per basis function
			q<-dim(object@U)[2]
			voxSel<-object@total2maskIndex[voxelIndex3D[1],voxelIndex3D[2],voxelIndex3D[3],1,1,1,1]
			
			#Get gamma coefficients of requested voxel and HRF basis function (which.beta)
			if(object@method=="bayes"){
				if(voxSel==0)return("Error: selected voxel is not within analysing mask.")
				coef<-apply(object@methodResults$gammaSamples[voxSel,,,which.beta],1,mean)
			}else{
				#Get element indices of betaHat
				if(which.beta>object@nStimBF){
					print("nStimBF < which.beta => take elements corresponding to 1st HRF BF")
					which.beta<-1
				}
				start<-1
				end<-p
				if(object@nStimBF>1 & which.beta==2){
					start<-p+1
					end<-2*p
				}
				if(object@nStimBF>2 & which.beta==3){
					start<-2*p+1
					end<-3*p
				}
				start<-start+q
				end<-end+q

				coef<-object@methodResults$thetaHat[voxelIndex3D[1],voxelIndex3D[2],voxelIndex3D[3],start:end,1,1,1]
			}
			#Calculate point estimate
			betaHatEst<-beta.hat.func(object@VC, VC.plot, coef,object@nVCknots,object@degreeVCspline)
			
			#Calculate interval estimate
			if(object@method=="bayes"){
				betaHatCi<-beta.hat.ci(object@methodResults$gammaSamples[voxSel,,,which.beta],object@nVCknots,object@degreeVCspline, object@VC, VC.plot, which.beta, alpha)
			}else{
				#calculate total design matrix
				Xtotal<-cbind(object@U,object@X)
				
				#Calculate total penalty matrix
				nCol=dim(Xtotal)[2]
				unpen<-diag(rep(1,q))
				pen<-diag(rep(1,object@nStimBF))%x%object@penMat
				penMatTotal<-matrix(0,nrow=nCol,ncol=nCol)
				penMatTotal[1:q,1:q]<-unpen
				penMatTotal[(q+1):nCol,(q+1):nCol]<-pen
				
				#Calculate CIs
				betaHatCi<-beta.penkq.ci(object@VC, VC.plot, betaHatEst,Xtotal, penMatTotal,object@nVCknots,object@degreeVCspline, object@methodResults$sigma2Hat[voxelIndex3D[1],voxelIndex3D[2],voxelIndex3D[3],1,1,1,1], which.beta, object@nStimBF, alpha)				
			}

			betaHat<-data.frame(cbind(betaHatEst,betaHatCi))
			names(betaHat)<-c("Beta","Lower","Upper")
			return(betaHat)
		}
)

#############################################################################
# Function to calculate point and interval estimates of beta_i(VC)
# for all voxels i=1,...,N -> and save them in a 4D nifti
#############################################################################

setGeneric("betaHatNifti",function(object,value,...){standardGeneric("betaHatNifti")})
setMethod(
		f="betaHatNifti",
		signature="ResultsVC",
		function(object,VC.plot=object@VC.plot,which.beta=1){
			dims<-dim(object@total2maskIndex)
			#Initialize nifti for all voxelwise betaHat-Estimates 
			betaHatTotal<-nifti.image.new()
			betaHatTotal$dim<-c(dims[1:3],length(object@VC.plot)) 
			nifti.image.alloc.data(betaHatTotal)
			
			for(i in 1:dims[1]){
				for(j in 1:dims[2]){
					for(k in 1:dims[3]){
						if(object@total2maskIndex[i,j,k,1,1,1,1]!=0){
							betaHatTotal[i,j,k,,1,1,1,1]<-as.vector(betaHat(object,voxelIndex3D=c(i,j,k),VC.plot=VC.plot,which.beta=which.beta)[,1])
						}
					}
				}
			}
			
			return(betaHatTotal)
		}
)

#############################################################################
# Function to calculate point and interval estimates of hrf_i(VC,t) for a
# specified value VC and t for a single voxel i.
# 
# Arguments: 
#	- object		object of class ResultsVC
#	- VCVal			value of VC variable for which hrf_i should be calculated
#	- tVal			time value (peristimulus time) for which hrf_i should be calculated
#	- voxelIndex3D	3D vector with 3D index of requested voxel			
#	- which.beta	corresponding HRF basis function number (default: canonical HRF <=> which.beta=1)
#############################################################################


setGeneric("hrfHatPoint",function(object,value,...){standardGeneric("hrfHatPoint")})
setMethod(
		f="hrfHatPoint",
		signature="ResultsVC",
		function(object,VCVal=median(object@VC),tVal,voxelIndex3D,alpha=object@alpha){
			p<-object@nVCknots+object@degreeVCspline-1 #Number of gamma parameters per basis function
			q<-dim(object@U)[2] #Number of nuisance parameters (confounders, baseline)
			voxSel<-object@total2maskIndex[voxelIndex3D[1],voxelIndex3D[2],voxelIndex3D[3],1,1,1,1] #1D index of selected voxel
			if(voxSel==0)return("Error: selected voxel is not within analysing mask.")
			
			#Get matrix with gamma coefficients of requested voxel and HRF basis function (which.beta)
			if(object@method=="bayes"){
				coefMat<-matrix(apply(object@methodResults$gammaSamples[voxSel,,,],c(1,3),mean),nrow=p,ncol=object@nStimBF)
			}else{
				start<-q+1
				end<-dim(object@methodResults$thetaHat)[4]
				coefVec<-object@methodResults$thetaHat[voxelIndex3D[1],voxelIndex3D[2],voxelIndex3D[3],start:end,1,1,1]
				coefMat<-matrix(coefVec,nrow=p,ncol=object@nStimBF)
			}
			hrfHatEst<-HRF.J.t.hat(object@VC, JVal=VCVal,tVal=tVal, coefMat,object@nVCknots,object@degreeVCspline)
			
			
			#Calculate interval estimate
			if(object@method=="bayes"){
				hrfHatCi<-hrf.J.t.hat.bayes.ci(object@VC, VCVal, tVal, object@methodResults$gammaSamples[voxSel,,,],object@nVCknots,object@degreeVCspline,alpha)
			}else{
				#calculate total design matrix
				Xtotal<-cbind(object@U,object@X)
				
				#Calculate total penalty matrix
				nCol=dim(Xtotal)[2]
				unpen<-diag(rep(1,q))
				pen<-diag(rep(1,object@nStimBF))%x%object@penMat
				penMatTotal<-matrix(0,nrow=nCol,ncol=nCol)
				penMatTotal[1:q,1:q]<-unpen
				penMatTotal[(q+1):nCol,(q+1):nCol]<-pen
				
				#Calculate CIs
				hrfHatCi<-hrf.J.t.hat.penkq.ci(object@VC,VCVal,tVal, coefMat,object@nVCknots,object@degreeVCspline, Xtotal, penMatTotal,  object@methodResults$sigma2Hat[voxelIndex3D[1],voxelIndex3D[2],voxelIndex3D[3],1,1,1,1], alpha)			
			}
			return(c(hrfHatEst,hrfHatCi))
		}
)


#############################################################################
# Function to calculate point and interval estimates of hrf_i(VC,t) for a
# specified value VC and t for a single voxel i.
# 
# Arguments: 
#	- object		object of class ResultsVC
#	- VC			value or vector of VC variable values for which hrf_i should be calculated
#	- t				time value or vector (peristimulus time) for which hrf_i should be calculated
#	- voxelIndex3D	3D vector with 3D index of requested voxel			
#	- alpha			significance/credibility level for pointwise confidence intervals
#############################################################################

setGeneric("hrfHat",function(object,value,...){standardGeneric("hrfHat")})
setMethod(
		f="hrfHat",
		signature="ResultsVC",
		function(object,VC=median(object@VC),t=object@t.plot,voxelIndex3D,alpha=object@alpha){
			if(is.matrix(VC))if(dim(VC)[2]>1)return("VC must be either a vector or a matrix with 1 column.")
			if(is.matrix(t))if(dim(t)[2]>1)return("t must be either a vector or a matrix with 1 column.")
			#Convert to matrix format 
			if(is.vector(VC))VCMat<-matrix(VC,nrow=length(VC),ncol=1)
			if(is.vector(t))tMat<-matrix(t,nrow=length(t),ncol=1)
			
			VCLen<-dim(VCMat)[1]
			tLen<-dim(tMat)[1]
			if(VCLen==1 & tLen==1){
				resMat<-hrfHatPoint(tVal=tMat,object=object,VCVal=VCMat,voxelIndex3D=voxelIndex3D,alpha=alpha)
				resMat<-matrix(resMat,ncol=3)
				colnames(resMat)<-c("HRF","Lower","Upper")
			}
			if(VCLen==1 & tLen>1){
				fcttmpOverT<-function(X){
					hrfHatPoint(tVal=X,object=object,VCVal=VCMat,voxelIndex3D=voxelIndex3D,alpha=alpha)
				}
				resMat<-as.data.frame(t(apply(tMat,1,fcttmpOverT)))
				colnames(resMat)<-c("HRF","Lower","Upper")
			}
			if(VCLen>1 & tLen==1){
				fcttmpOverVC<-function(X){
					hrfHatPoint(tVal=tMat,object=object,VCVal=X,voxelIndex3D=voxelIndex3D,alpha=alpha)
				}
				resMat<-as.data.frame(t(apply(VCMat,1,fcttmpOverVC)))
				colnames(resMat)<-c("HRF","Lower","Upper")
			}
			if(VCLen>1 & tLen>1){
				fcttmpOverVCt<-function(X){
					VCVal=X[1]
					tVal=X[2]
					hrfHatPoint(tVal=tVal,object=object,VCVal=VCVal,voxelIndex3D=voxelIndex3D,alpha=alpha)
				}				
				VCtgrid<-expand.grid(x=VCMat,t=tMat)
				resMatTmp<-apply(VCtgrid,1,fcttmpOverVCt)
				resMat<-array(0,dim=c(VCLen,tLen,3))
				resMat[,,1]<-matrix(resMatTmp[1,],nrow=VCLen,ncol=tLen)
				resMat[,,2]<-matrix(resMatTmp[2,],nrow=VCLen,ncol=tLen)
				resMat[,,3]<-matrix(resMatTmp[3,],nrow=VCLen,ncol=tLen)
			}
			
			return(resMat)
		}
)


########################################################################################################
#Show
setMethod(
		f="show",
		signature="ResultsVC",
		function(object){
			cat("This is an object of Class 'ResultsVC'. Package 'RfmriVC'.\n")
			cat("It contains the results and needed configuration information of a call to the RfmriVC-algorithm.\n\n")
			
			cat(paste("method = Requested algorithm:",object@method,"\n"))
			cat(paste("nStimBF = Number of canonical HRF Basis functions:",object@nStimBF,"\n"))
			cat(paste("For each effect of the ",object@nStimBF," HRF basis functions a VC was fitted.\n\n"))
			
			cat(paste("VC = Variable the varying coefficient depends on : "))
			cat(paste(head(object@VC)))
			cat(" ...\n")
			cat(paste("nVCknots = Number of knots for spline basis of VC:",object@nVCknots,"\n"))
			cat(paste("degreeVCspline = Degree of the VC spline:",object@degreeVCspline,"\n"))
			
			cat(paste("stimTimes = Vector of stimulus presentation times: "))
			cat(paste(head(object@stimTimes)))
			cat(" ...\n")

			
			
			cat(paste("VC.plot = Default VC vector for plotting: "))
			cat(paste(head(object@VC.plot)))
			cat(" ...\n")
			
			cat(paste("t.plot = Default t vector for plotting: "))
			cat(paste(head(object@t.plot)))
			cat(" ...\n")

			cat(paste("alpha = Default significance/credibility level:",object@alpha,"\n"))
			
			cat("Further slots: U (nuisance regressors), X (VC design matrix part), penMat (penalty matrix for one VC spline), \n total2maskIndex (mapping from total to mask index)\n")
			
			cat("Results can be found in slot 'methodResults'\n")
			
		}
)
