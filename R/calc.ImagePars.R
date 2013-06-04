calc.ImagePars <-
function(dims,zVec=NULL,zRange=NULL,mfrowVal=NULL){
  #Set zVec, i.e. vector of indices of image layers to be plotted
    #---------------------------------------------------------------
    #Distinction of cases:
    # -1- if neither zRange nor zVec is specified: Select all layers 
    if(is.null(zRange) & is.null(zVec)){
        zVec<-seq(1,dims[3])
        nlayers<-dims[3]
    }
    
    # -2- if zRange is specified (not zVec): Select all layers between zRange[1] and zRange[2]
    if(!is.null(zRange) & is.null(zVec)){
      #Set range of layers to be plotted
      if(zRange[1]<=0)zRange[1]<-1
      if(zRange[2]==0 | zRange[2]>dims[3])zRange[2]<-dims[3]
      zVec<-seq(zRange[1],zRange[2])
      nlayers<-zRange[2]-(zRange[1]-1)
    }
    
    # -3- if zVec is specified (regardless of zRange being specified): Selection of layers is already given in zVec
    # Check if selection is ok
    if(!is.null(zVec)){ #is.null(zRange) & 
      zlen<-length(zVec)
      if(zVec[1]<=0 | zVec[1]>dims[3]) zVec<-1:zVec[zlen]
      if(zVec[zlen]<=0 | zVec[zlen]>dims[3]) zVec<-zVec[1]:dims[3]
      nlayers<-zlen#zVec[zlen]-(zVec[zlen]-1)
    }
		

	
    # If not specified, calculate mfrow-dimensions for plotting
    #------------------------------------------------------------
	if(is.null(mfrowVal)){
		nRow<-0
		nCol<-0
		if(nlayers>4){
			nRow<-3
			nCol<-ceiling(nlayers/nRow)
			while((nRow/nCol)<0.55){
				nRow<-nRow+1
				nCol<-ceiling(nlayers/nRow)        
			}
			if(nRow>nCol){ #Prefer lanscape -> switch nRow and nCol
				ztmp<-nCol
				nCol<-nRow
				nRow<-ztmp
			}
		}else{
			nCol<-nlayers
			nRow<-1
		}
	}else{
		nRow<-mfrowVal[1]
		nCol<-mfrowVal[2]
		if(nRow*nCol<=nlayers){
			nRow<-ceiling(nlayers/nCol) 
		}
		
	}
    
    return(list(zVec=zVec,nRow=nRow,nCol=nCol))
}
