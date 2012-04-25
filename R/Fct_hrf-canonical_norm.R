################################################################################
# Author: Stefanie Kalus, Ludwig Bothmann
# Last change: 19.04.12
################################################################################



##############################################################################
# Basis functions for canonical HRF set - in analogy to SPM implementation   #
##############################################################################

#Canonical HRF:
#---------------
cHRF<-function(x,dtp=0.125,p1=6,p2=16,p3=1,p4=1,p5=6,p6=0){
  y <- dgamma(x/dtp-p6,shape=p1/p3,scale=p3/dtp)-dgamma(x/dtp-p6,shape=p2/p4,scale=p4/dtp)/p5 
  y <- y/0.02193015  #scale to (Max(abs)) = 1
  return(y)
}

#Time derivative:
#------------------
# Derivative component
dgammadp6<-function(x,dtp=0.125,p1=6,p3=1,p6=0){
  y<-rep(NA,length(x))
  for(i in 1:length(x)){
    if(x[i]<=0){ y[i]<-0}
    else{y[i]<-(dtp/p3)^(p1/p3)/gamma(p1/p3)*(x[i]/dtp-p6)^((p1/p3)-2)*exp(-dtp/p3*(x[i]/dtp-p6))*(-(p1/p3-1)+ (x[i]/dtp-p6)*(dtp/p3))}
  }  
  return(y)
}

#Whole time derivative
cHRFdp6<-function(x,dtp=0.125,p1=6,p2=16,p3=1,p4=1,p5=6,p6=0){
  	y<-dgammadp6(x,dtp,p1,p3,p6)-dgammadp6(x,dtp,p2,p4,p6)/p5
	y<-y/0.001071197		#scale to (Max(abs)) = 1
}

#Dispersion derivative:
#-----------------------
cHRFdp3<-function(x,dtp=0.125,p1=6,p3=1,p6=0){
  y<-rep(NA,length(x))
  for(i in 1:length(x)){
    if(x[i]<=0){ y[i]<-0}
    else{ y[i]<- - (dtp/p3)^(p1/p3)/gamma(p1/p3)*(x[i]/dtp-p6)^(p1/p3-1)*exp(-dtp*(x[i]/dtp-p6)/p3)*(p1*log(dtp/p3)+p1-digamma(p1/p3)*p1+p1*log(x[i]/dtp-p6)-dtp*(x[i]/dtp-p6))/p3^2}
  }
  
  y <- y/0.01061682  #scale to (Max(abs)) = 1

  return(y)
}

