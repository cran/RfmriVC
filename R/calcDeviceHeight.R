calcDeviceHeight <-
function(widthVal,mfrowVal,maiVal,omiVal,niftiDim){
   nRow<-mfrowVal[1]
   nCol<-mfrowVal[2]
   dx<-(widthVal-(omiVal[2]+omiVal[4]))/nCol- (maiVal[2]+maiVal[4])  #width of one figure
   dy<-niftiDim[2]/niftiDim[1]*dx  #Get height of one figure according to dim-ratio
   heightVal<-(maiVal[1]+maiVal[3] + dy)*nRow + omiVal[1]+omiVal[3]
   return(heightVal)
}
