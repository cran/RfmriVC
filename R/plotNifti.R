plotNifti <-
function(nifti.img,zRange=NULL,zVec=NULL,mfrowVal=NULL,plotTitle=TRUE,pdfoutput=FALSE,pdfName=NULL,scale=1,...){
    require(Rniftilib)
    dims<-dim(nifti.img)
    
    #Calculate zVec with slice indices to be plottet
    # and nRow and nCol 
    #-----------------------
    imagePars<-calc.ImagePars(dims,zRange=zRange,zVec=zVec,mfrowVal)
    zVec<-imagePars$zVec
    nRow<-imagePars$nRow
    nCol<-imagePars$nCol
	mfrowVal<-c(nRow,nCol)

    #Extract image name
    #----------------------
    if(nifti.img$fname!=""){
      strparts1a<-strsplit(nifti.img$fname,"\\",fixed=TRUE)[[1]]
      strparts1b<-strsplit(strparts1a[length(strparts1a)],".",fixed=TRUE)[[1]]
	  #titleString<-c(paste(strparts1b[1]))
      titleString<-c(paste(strparts1b[1],"\n layers:",paste(zVec,collapse=" ")))
    }else{titleString<-""}
        
    #If requested: open pdf
    #------------------------
    
    #Set pars:
	tMa<-as.numeric(plotTitle)*0.65 #title Margin
	maiVec<-c(0,0, 0,0)+0.015
	omiVec<-c(0,0,tMa,0)+0.15
    cex.main.val<-1
    
    #Calculate appropriate width and height, so that image ratios are preserved
    widthVal<-7
	heightVal<-calcDeviceHeight(widthVal,mfrowVal,maiVec,omiVec,dims)
	#print(heightVal)
    if(pdfoutput){
		if(is.null(pdfName)){
			pdf(paste(strparts1b[1],".pdf",sep=""),width=widthVal,height=heightVal)
		}else{
			pdf(paste(pdfName,".pdf",sep=""),width=widthVal,height=heightVal)
		}
	}
	
    par(mfrow=mfrowVal,mai=maiVec,omi=omiVec,...)
    
    zMin<-min(nifti.img[,,,1,1,1,1]/scale)
    zMax<-max(nifti.img[,,,1,1,1,1]/scale)
    colormap<- gray(1:1000/1000)
	
    #Plot image layers
    for(i in zVec){
      image(nifti.img[,,i,1,1,1,1]/scale,zlim=c(zMin,zMax),col=colormap,axes=FALSE) #ylab=paste("z=",i)
      #box()
    }
    #if(plotTitle)title(nifti.img$fname,outer=T)
	if(plotTitle)title(titleString,outer=T)
    if(pdfoutput)dev.off()
}
