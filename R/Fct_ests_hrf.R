################################################################################
# Point and interval estimates for voxelwise hrf estimates
# Created on: 15.11.11
# Last change: 24.04.12
# Author: Ludwig Bothmann
# Changed by: Stefanie Kalus
################################################################################

#############################################################################
# Mit der Funktion HRF.hat() soll die geschätzte HRF berechnet werden.
# Übergeben werden muss:
# 	- J.vig: 	Zeitreihe des Vigilanzparameters (Zu jedem Stimuluszeitpunkt ein Wert)
#	- JVal:		Der Wert von J, für den die HRF.hat berechnet werden soll
# 	- tVal:		t-Wert, für die die HRF berechnet werden soll
# 	- coefMat: 	Matrix (nGamma x nStimBF) der geschätzten Gamma Parameter
# 	- K: 		Anzahl der Knoten für die Spline-Basis
# 	- l: 		Grad der Spline-Basis
#############################################################################

HRF.J.t.hat <- function(J.vig, JVal=median(J.vig),tVal=0, coefMat, K, l){
	nStimBF<-dim(coefMat)[2]
	
	# Berechnung von beta.hat
	beta.hat <- c(beta.hat.func(J.vig,JVal,coefMat[,1],K,l),0,0)
	if(nStimBF>1){
		beta.hat[2]<-beta.hat.func(J.vig,JVal,coefMat[,2],K,l)
	}
	if(nStimBF>2){
		beta.hat[3]<-beta.hat.func(J.vig,JVal,coefMat[,3],K,l)
	}
	
	# Werte der 3 kanonischen Basisfunktionen zum t-Grid
	b1 <- cHRF(tVal)
	b2 <- cHRFdp6(tVal)
	b3 <- cHRFdp3(tVal)
	
	bVec<-c(b1,b2,b3)
	
	# Die geschätzte HRF ist eine gewichtete Summe, siehe 1.2.1 auf Seite 3
	# hrf.hat ist ein Vektor, der die gleiche LÃ¤nge hat wie t.g
	hrf.hat <- sum(bVec*beta.hat)
	
	return(hrf.hat)
}


################################################################################
#								BAYES
# Mit der Funktion hrf.J.t.hat.bayes.ci() werden die Grenzen für ein 
# Kredibilitätsintervall für die HRF zu Zeitpunkt t (PST) und J-Wert JVal 
# im Modell mit variierendem Stiumuluseffekt bestimmt.
#
# Übergeben werden muss:
# 	- J.vig: 	Zeitreihe des Vigilanzparameters (Zu jedem Stimuluszeitpunkt ein Wert)
#	- JVal:		Der Wert von J, für den die HRF.hat berechnet werden soll 
# 	- tVal:		t-Wert, für die die HRF berechnet werden soll
#	- samples: 	Eine Matrix mit 3*p Spalten, welche die Samplingpfade der 
#						gammas enthalten
# 	- alpha:   	Ein Vektor der angibt, welche Grenzen von Interesse sind
#
################################################################################

hrf.J.t.hat.bayes.ci <- function(J.vig, JVal, tVal, samples, K, l, alpha, ...){
	ifelse(is.na(dim(samples)[3]),nStimBF<-1,nStimBF<-dim(samples)[3])
	
	# Bestimmen der 3 Basisfunktionen
	tVal<-as.vector(tVal)
	B1 <- cHRF(tVal)
	B2 <- cHRFdp6(tVal)
	B3 <- cHRFdp3(tVal)
		
	# Bestimmen der Werte der B-Spline-Basisfunktionen zur gegebenen Vigilanzzeitreihe
	# an der Stelle JVal
	e.J <- matrix(E.design(J.vig, K, l, J.plot=JVal),nrow=1)
	
	# Initialisieren der Matrix, in die die Samples der beta_k(J) geschrieben werden

	# Bestimmen der Samples der beta_k(J)
	ifelse(nStimBF==1, samples.hrf<-B1*e.J %*% samples , samples.hrf <-B1*e.J %*% samples[,,1])
	if(nStimBF>1) samples.hrf <- rbind(samples.hrf,B2*e.J %*% samples[,,2])
	if(nStimBF>2) samples.hrf <- rbind(samples.hrf,B3*e.J %*% samples[,,3])
	samples.hrf<-apply(samples.hrf,2,sum)
	quants<-quantile(samples.hrf,probs=c(alpha/2,1-alpha/2))
		
	return(quants)	
}


################################################################################
#								Classical
# Mit der Funktion hrf.J.t.hat.bayes.ki() werden die Grenzen für ein 
# Kredibilitätsintervall für die HRF zu Zeitpunkt t (PST) und J-Wert JVal 
# im Modell mit variierendem Stiumuluseffekt bestimmt.
#
# Übergeben werden muss:
# 	- J.vig: 	Zeitreihe des Vigilanzparameters (Zu jedem Stimuluszeitpunkt ein Wert)
#	- JVal:		Der Wert von J, für den die HRF.hat berechnet werden soll
# 	- tVal:		t-Wert, für die die HRF berechnet werden soll
# 	- coefMat: 	Matrix (nGamma x nStimBF) der geschätzten Gamma Parameter
# 	- K: 		Anzahl der Knoten für die Spline-Basis
# 	- l: 		Grad der Spline-Basis
#	- X:			Die gesamte Designmatrix
#	- Pen:			Die Penalisierungsmatrix für's gesamte Modell
#	- sigma2.hat:	Die geschätzte Residuenvarianz aus dem LM
# 	- alpha:		Das Signifikanzniveau
#	
################################################################################


hrf.J.t.hat.penkq.ci <- function(J.vig,JVal,tVal, coefMat, K, l, X, Pen,  sigma2.hat, alpha){
	nStimBF<-dim(coefMat)[2]
	hrf<-HRF.J.t.hat(J.vig, JVal=JVal,tVal=tVal, coefMat, K, l)
	
	# Bestimmen der 3 Basisfunktionen
	B1 <- cHRF(tVal)
	B2 <- cHRFdp6(tVal)
	B3 <- cHRFdp3(tVal)
	
	# Bestimmen des Matrix B
	B <- matrix(c(B1,B2,B3),nrow=1,ncol=3)
	B <- matrix(B[,1:nStimBF],nrow=1,ncol=nStimBF)
	
	
	# Bestimmen der Kovarianzmatrix von beta.dach.pen.total
	Kov <- sigma2.hat * solve(t(X) %*% X + Pen) %*% t(X) %*% X %*% solve(t(X) %*% X + Pen)
	
	# Anzahl Basisfunktionen
	p <- K+l-1
	
	# Von dieser Kovarianzmatrix brauche ich nur Block für gamma
	# Gesamter gamma-Block
	Kov.gamma <- Kov[(nrow(Kov)-nStimBF*p + 1):nrow(Kov),(nrow(Kov)-nStimBF*p + 1):nrow(Kov)]
	
	# Bestimmen von E.stern
	E.stern <- E.stern.design(J.vig, K, l, JVal,nStimBF)
	
	# Bestimmen der Varianz der HRF
	var.hrf<- B %*% E.stern %*% Kov.gamma %*% t(E.stern) %*% t(B)
	
	#Ci
	ci<- hrf + c(-1,1)*qnorm(1-alpha/2)*sqrt(var.hrf)
	
	
	return(ci)
}

