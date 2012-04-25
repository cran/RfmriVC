################################################################################
# Functions to calculate design matrix parts of fMRI VC model
# Created on: 14.11.11
# Last change: 19.04.12
# Author: Ludwig Bothmann
# Changed by: Stefanie Kalus
################################################################################




######################################################################
# HPF basis functions (DCT set)
######################################################################

# 1st basis function equals intercept
#-------------------------------------
w1 <- function(t){
	
	# Bestimme die Länge von t
	T <- length(t)
	
	# Initialisiere Vektor, der später zurückgegeben wird
	w <- vector(length=T)
	
		# Befülle den Vektor mit der Konstanten 1/Wurzel(2)
		for(i in 1:T){
			w[i] <- 1/(sqrt(2))
		}
	
	return(w)
}

# k-th basis function is cosine/sine oscillation
#------------------------------------------------
wk <- function(k, t) {
	
	# Bestimme die Länge von t
	T <- length(t)

	# Initialisiere Vektor, der später zurückgegeben wird
	w <- vector(length=T)
	
		# Befülle den Vektor mit den ensprechenden Werten gemäÃŸ der Formel
		for(i in 1:T){
			w[i] <- sqrt(2/T) * cos((pi*(2*t[i]-1)*(k-1))/(2*T))
		}
	
	return(w)
}


# Calculate whole HPF design matrix W = cbind(HPF basis functions)
#------------------------------------------------------------------

w.design <- function(t.grid,freqHPF,tr){
	
	T <- length(t.grid)
	
	# Berechnung der Anzahl Basisfunktionen
	p.1 <- floor(4*T / (2*freqHPF/tr - 1) + 1) # => 10
	
	# if(T != 354){
	# 	print("Parameter p_1 der DCT-Baseline überprüfen")
	# }
	
	
	W <- matrix(nrow=T,ncol=p.1)
	W[,1] <- w1(t.grid)


	for(i in 2:p.1){
		W[,i] <- wk(i,t.grid)
	}

	W <- data.frame(W)
	
	return(W)
}



###########################################################################
# Vektor B_t für den Zeitpunkt t einer fMRT-Messung berechnen. Übergeben werden muss:
#	- stimTimes: 	Zeitreihe der Stimuluszeitpunkt
#	- t: 	Index der interessierenden fMRT-Messung
#	- tr:	Repetitiontime = Zeit zwischen zwei fMRT-Messungen, default=2
###########################################################################

B.design <- function(stimTimes, t, tr, nStimBF,...){
		
	# Anzahl Stimuli
	M <- length(stimTimes)
	
	# Realzeit der fMRT-Messung
	time <- tr*(t-1)
	
	# Vektor B initialisieren
	B <- vector(length= nStimBF*M)
	
	# 1. Basisfunktion: Werte zur Zeitdifferenz von fMRT-Messung und Stimulus 
		# pro Stimulus berechnen
	for(i in 1:M){
		B[i] <- cHRF(time - stimTimes[i],...)
	}
	
	# 2. Basisfunktion analog
	if(nStimBF>1){
		for(i in 1:M){
			B[ (M+i) ] <- cHRFdp6(time - stimTimes[i],...)
		}
	}

	# 3. Basisfunktion analog
	if(nStimBF>2){
		for(i in 1:M){
			B[ (2*M+i) ] <- cHRFdp3(time - stimTimes[i],...)
		}
	}
	
	return(B)
	
}

######################################################################
# Matrix E.stern(J) berechnen. Übergeben werden muss:
#	- J.vig:	Zeitreihe des Vigilanzparameters (Zu jeder fMRT-Messung ein Wert)
# 	- K: 		Anzahl der Knoten für die Spline-Basis
# 	- l: 		Grad der Spline-Basis		
#	- nStimBF:	Anzahl an kanonischen HRF Basisfunktionen
#  	- optional: 	Wert J.plot für den gezeichnet werden soll. 
#		Diese Erweiterung wurde nötig, da ich jetzt quantilsbasierte Knoten wähle. 
#		Wenn man nicht plotten will, sondern die Designmatrix bestimmen, 
#		dann kann dieser Parameter unspezifiziert bleiben. (Beim Plotten
#		brauche ich J.vig also für die Bestimmung der Knotenpunkte, bei 
#		der Bestimmung der Designmatrix auch und zusätzlich um die Aufmerksamkeit
#		beim Stimuluszeitpunkt zu bekommen.)
######################################################################

E.stern.design <- function(J.vig, K, l, J.plot=J.vig,nStimBF){
	
	M <- length(J.plot)
	
	# Anzahl zu berechnender Parameter
	p <- K+l-1	
	
	E <- E.design(J.vig, K, l, J.plot)
	if(nStimBF==1){
		E.stern<-E
	}else{
		E.stern <-diag(rep(1,nStimBF))%x%E  
	}
	return(E.stern)
}

######################################################################
# Matrix E(J) berechnen. Diese Funktion berechnet die Matrix der B-Spline-Basis.
#	Es ergeben sich M Zeilen und p Spalten.
#	Übergeben werden muss:
#	- J.vig:	Zeitreihe des Vigilanzparameters (Zu jeder fMRT-Messung ein Wert)
#  	- K:		Anzahl der Knoten
#  	- l:		Grad der Spline-Basis
#	- optional: J.plot
######################################################################



E.design <- function(J.vig, K, l, J.plot=J.vig){
			
	knots <- quantile(J.vig, seq(1,K-2,by=1) / (K-1) )
	
	for(i in 1:length(J.plot))if((J.plot[i]>max(J.vig)) | (J.plot[i] <min(J.vig))){print(paste(J.plot,"not in range J.vig: (",min(J.vig),"),(",max(J.vig),")"))}

	E<-bs(x=J.plot, knots=knots, degree=l, intercept=TRUE,Boundary.knots=range(J.vig)) # ist intercept=TRUE eine gute Wahl?
		# Update 21.12.11: Vielleicht ist doch FALSE besser wegen Identifizierbarkeit, ich probier es mal...
		# Update: 22.12.11: Intercept=TRUE müsste schon passen, sonst fängt beta immer bei 0 an
	

	return(E)
}



######################################################################
# Designmatrix für den Teil der Stimuluseffekte: X berechnen.
#	Übergeben werden muss: 
#	- T:		Anzahl an fMRT-Messungen
#	- stimTimes: 		Zeitreihe der Stimuluszeitpunkt
#	- tr:		Repetitiontime = Zeit zwischen zwei fMRT-Messungen, default=2
#	- J.vig:	Zeitreihe des Vigilanzparameters (Zu jeder fMRT-Messung ein Wert)
#  	- K:		Anzahl der Knoten
#  	- l:		Grad der Spline-Basis
######################################################################

x.design <- function(T, stimTimes, tr, J.vig, K, l,nStimBF){

	# Anzahl zu berechnender Parameter
	p <- K+l-1	
		
	# Matrix X initialisieren
	X <- matrix(0,nrow=T, ncol=(nStimBF*p))
	
	# Matrix E.stern berechnen
	E.stern <- E.stern.design(J.vig, K, l,J.vig,nStimBF)
	
	
	for(t in 1:T){
	
		# Vektor B berechnen
		B_t <- B.design(stimTimes, t, tr, nStimBF)
		
		# Zeilenweise befüllen mit B_t^T %*% E.stern(J)
		X[t,] <- t(B_t) %*% E.stern
		
	}
	
	return(X)
}




######################################################################
# Penaltymatrix
######################################################################

k.pen <- function(n,grad){
	
	a <- diag(1,n-1)
	b <- diag(-1,n-1)
	a1 <- cbind(rep(0,n-1),a)
	b1 <- cbind(b,rep(0,n-1))
	d <- a1+b1
	t(d)%*%d
	d2 <- d[1:(n-2),1:(n-1)]%*%d
	d3 <- d2[1:(n-3),1:(n-1)]%*%d
	k1 <- t(d)%*%d
	k2 <- t(d2)%*%d2
	k3 <- t(d3)%*%d3

	ifelse(grad==1,return(k1),ifelse(grad==2,return(k2),return(k3)))
	
}










