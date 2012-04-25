################################################################################
# Point and interval estimates for spatially varying beta coefficients
# Created on: 15.11.11
# Last change: 24.04.12
# Author: Ludwig Bothmann
# Changed by: Stefanie Kalus
################################################################################

#######################################################################################
# Mit der Funktion beta.hat.func() sollen die Schätzer beta.hat, d.h. Schätzungen für 
#	beta_k(J)=E(J)%*%gamma_k 
# für HRF Basisfunktion k in Abhängigkeit von J berechnet werden.
# Übergeben werden muss: 
# 	- J.vig: 	Zeitreihe der beobachteten Variieredned Koeff. Variable J (=VC) (Zu jedem Stimuluszeitpunkt ein Wert)
# 	- J.plot: 	Vektor mit J-Werten, für die beta.hat berechnet werden soll
# 	- coef: 	Vektor der geschätzten gamma_k Parameter für HRF Basisfunktion k
# 	- K: 		Anzahl der Knoten für die Spline-Basis
# 	- l: 		Grad der Spline-Basis
########################################################################################

beta.hat.func <- function(J.vig, J.plot, coef, K, l){
	
	# Matrix der Spline-Basis	
	e.stern <- E.stern.design(J.vig, K, l, J.plot,nStimBF=1)
	
	# Vektor der geschätzten Koeffizienten
	gamma.hat <- coef
	
	beta.hat <- e.stern%*%gamma.hat
	
	return(beta.hat)	
}

################################################################################
#								BAYES
# Mit der Funktion beta.hat.ci() werden die Grenzen für ein gleichendiges
# Kredibilitätsintervall für ein beta_k(J) im Modell mit variierendem 
# Stiumuluseffekt bestimmt. 
#
# Übergeben werden muss:
#	- samples:		Eine Matrix, die die Samplingpfade von Parametervektor gamma_k enthält
#	- K:			Die Anzahl der Knoten für die B-Spline-Basis
#	- l:			Der Grad der Splines
#	- J.vig:		Die Vigilanzzeitreihe
#	- J.plot:		Ein Vektor, der angibt, für welche J die Grenzen bestimmt werden sollen
#	- which.beta: 	Index der angibt, für welches beta_k die Grenzen bestimmt werden sollen
# 	- alpha:		Ein Vektor der angibt, welche Grenzen von Interesse sind
#	
# Zurückgegeben wird eine Matrix, die in den Zeilen die gewünschten 
#	Quantile für den jeweiligen Wert von J.plot enthält.
################################################################################


beta.hat.ci <- function(samples, K, l, J.vig, J.plot, which.beta, alpha){
		
	# Bestimmen der Werte der B-Spline-Basisfunktionen zur gegebenen Vigilanzzeitreihe
	# an der Stelle J.plot
	e.J <- E.design(J.vig, K, l, J.plot)
		
	# Bestimmen der Samples des interessierenden beta_k(J)
	samples.beta <- e.J%*%samples
	
	# Bestimme alpha Quantile
	quantsLower <-apply(samples.beta,1,quantile,probs=alpha/2)
	quantsUpper <-apply(samples.beta,1,quantile,probs=1-alpha/2)
	quants<-cbind(matrix(quantsLower,ncol=1),matrix(quantsUpper,ncol=1))
	
	return(quants)
	
}





################################################################################
#								CLASSICAL
# Mit der Funktion beta.penkq.ci() werden die Grenzen für ein gleichendiges
# Konfidenzintervall für die HRF im Modell mit variierendem 
# Stimuluseffekt bestimmt. 
#
# Übergeben werden muss:
#	- J.vig:	   	Vektor der Vigilanzzeitreihe J(t)
# 	- J.plot:  		Der Wert von J, für den die KIs der HRF bestimmt werden sollen
#	- beta:			Der Vektor mit den Werten des geschätzten beta_k-Vektors
#	- X:			Die (gesamte) Designmatrix
#	- Pen:			Der Penalisierungsterm
#	- K:			Die Anzahl der Knoten für die B-Spline-Basis
#	- l:			Der Grad der Splines
#	- sigma2.hat:	Die geschätzte Residuenvarianz aus dem LM
#	- which.beta: 	Index der angibt, für welches beta_k die Grenzen bestimmt werden sollen
# 	- alpha:		Das Signifikanzniveau
#	
#
# Zurückgegeben wird eine Matrix, die in den Zeilen die gewünschten 
#	Quantile für den jeweiligen Wert von J.plot enthält.
################################################################################

beta.penkq.ci <- function(J.vig, J.plot, beta, X, Pen, K, l, sigma2.hat, which.beta, nStimBF, alpha){
	

	# Bestimmen der Kovarianzmatrix von beta.dach.pen.total
	Kov <- sigma2.hat * solve(t(X) %*% X + Pen) %*% t(X) %*% X %*% solve(t(X) %*% X + Pen)
	
	# Anzahl Basisfunktionen
	p <- K+l-1
	
	# Von dieser Kovarianzmatrix brauche ich nur Block für gamma
	# Gesamter gamma-Block
	Kov.gamma <- Kov[(nrow(Kov)-nStimBF*p + 1):nrow(Kov),(nrow(Kov)-nStimBF*p + 1):nrow(Kov)]
	
	
	
	# Nur Block von gamma_k
	# <-> entsprechend which.beta und nStimbf
	if(which.beta>nStimBF)print("nStimBF < which.beta => take elements corresponding to 1st HRF BF")
	start<-1
	end<-p
	if(nStimBF>1 & which.beta==2){
		start<-p+1
		end<-2*p
	}
	if(nStimBF>2 & which.beta==3){
		start<-2*p+1
		end<-3*p
	}
	Kov.gamma_k <- Kov.gamma[start:end,start:end]

	
	# Bestimmen der Werte der B-Spline-Basisfunktionen zur gegebenen Vigilanzzeitreihe
	# an den Stellen von J.plot
	e.J <- t(E.design(J.vig, K, l, J.plot))	
	
	# Bestimmen der Varianz der HRF
	var.beta <- vector(length=length(J.plot))
	
	for(i in 1:length(J.plot)){
		var.beta[i] <- t(e.J[,i]) %*% Kov.gamma_k %*% e.J[,i]
	}
	
	
	# Matrix der KIs initialisieren 
	ci <- matrix(nrow=length(J.plot), ncol=2)

	# Matrix der KIs mit den Werten befüllen
	for(i in 1:length(J.plot)){
		# Bestimmen des KIs
		ci[i,] <- beta[i] + c(-1,1)*qnorm(1-alpha/2)*sqrt(var.beta[i])
	}

	return(ci)
}


#############################################################################
# Function to plot estimate of beta(J) against beta with KIs
# Arguments: 
# 	- beta.hat.est	matrix with 1st column betaHat 2nd and 3rd column upper and lower KI
#	- main:		(optional) non-default plot title
#############################################################################

beta.hat.plot <- function(VC.plot,beta.hat.est,xlab="VC",ylab=expression(beta~(VC)), main="",...){
	ylimVec<-c(min(beta.hat.est),max(beta.hat.est))
	par(mar=c(5, 6, 4, 2) + 0.1,...)
	plot(VC.plot,beta.hat.est[,1],type="l",xlab=xlab,ylab=ylab,ylim=ylimVec,main=main)
	lines(VC.plot,beta.hat.est[,2])
	lines(VC.plot,beta.hat.est[,3])
}
