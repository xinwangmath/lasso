# The lasso modification of the LAR algorithm

# input: myY -- response variable; myX -- predictors
# output: coefficients, predicted values for every step

source("lassoClass.R")


lasso = function(myY, myX) {
	n = length(myY)
	#print(class(myX))
	p = dim(myX)[2]
	predictorNames = attributes(myX)$names

	if(class(myX) != "data.frame"){
		myX = data.frame(myX)
	}


	# standardize myX 
	myXStar = apply(myX, 2, function(xx){
		return( (xx - mean(xx))/sd(xx) )
		})
	
	myXStar = data.frame(myXStar)
	mu = matrix(mean(myY), n, 1)

	currentList = c()
	potentialList = attributes(myXStar)$names

	betaHat = matrix(0, p, min(p, n-1))

	np = 0
	k = 1

	while( np < (min(p, n-1)) ){

		lassoFlag = FALSE

		if(k > min(p, n-1)){
			betaHat = cbind(betaHat, rep(0, p))
		}

		rp = myY - mu[, dim(mu)[2]]
		if(length(potentialList)>0){
			XPotential = myXStar[potentialList]
			j = which.max( abs( cor(XPotential, rp ) ) )
	

			oldCurrentList = currentList
			currentList = c(currentList, potentialList[j])
			potentialList = potentialList[-j]  
		}


		activeIndex = match(currentList, attributes(myXStar)$names)
		inactiveIndex = match(potentialList, attributes(myXStar)$names)
		oldActiveIndex = match(oldCurrentList, attributes(myXStar)$names)

		AkCorSign = as.vector(sign( t(as.matrix(myXStar[currentList])) %*% rp ) )
	

		designMat = as.matrix(myXStar[currentList]) %*% diag(AkCorSign, nrow = length(AkCorSign), ncol = length(AkCorSign))


		if(length(potentialList)>0){
			cHatkm = t(as.matrix(myXStar[potentialList])) %*% rp
			} else{
				cHatkm = c()
			}
		


		capitalCHat = unique( (t(designMat) %*% rp) )
		uniqueTor = 0.0001
		if(length(capitalCHat) > 1){
			if((max(capitalCHat) - min(capitalCHat)) > uniqueTor){
				print("errror, non unique capitalCHat")
			}
		}
		capitalCHat = capitalCHat[1]

		oneAk = rep(1, length(currentList))
		
		if(length(potentialList)>0){
			aVec = t(as.matrix(myXStar[potentialList])) %*% designMat %*% solve( (t(designMat) %*% designMat), oneAk )
			} else{
				aVec = c()
			}
		


		tempVec1 = (capitalCHat - cHatkm)/(1 - aVec)
		tempVec2 = (capitalCHat + cHatkm)/(1 + aVec)
		

		tempVecPositive = c(tempVec1[tempVec1>0], tempVec2[tempVec2>0])
		stepc = min(tempVecPositive)
		

		if( k > 1) {
			AkCorSign = as.vector(sign( t(as.matrix(myXStar[currentList])) %*% rp ) )
		    designMat = as.matrix(myXStar[currentList]) %*% diag(AkCorSign)

			bTilde = solve((t(designMat) %*% designMat), oneAk)
			dVec = diag(AkCorSign) %*% bTilde

			


			dVecTotal = rep(0, p)
			dVecTotal[activeIndex] = dVec

			
			newcVec = rep(-1, p)
			
			newcVec[activeIndex] = -betaHat[activeIndex, k-1]/dVec

			
			cStar = min(newcVec[newcVec > 0])
			
			

			if((cStar < stepc) && (cStar != Inf) ){
				stepc = cStar
				
				lassoFlag = TRUE

			}

		}

		# update the prediction muHat
		uTilde = stepc * designMat %*% solve( (t(designMat) %*% designMat), oneAk )
		muk = mu[, dim(mu)[2]] + uTilde
		mu = cbind(mu, muk)
		bVec = stepc * solve( (t(designMat) %*% designMat), oneAk )
		
		# this step is to restore the sign for betaHat
		bVec = AkCorSign * bVec
		


		oldInNewIndex = match(oldCurrentList, currentList)
		oldInNewLogical = oldCurrentList %in% currentList


		if(any(oldInNewLogical)){
			oldBetaHatIncrement = bVec[oldInNewIndex]
			newBetaHatIncrement = bVec[-oldInNewIndex]
		} else{
				newBetaHatIncrement = bVec
		}

		
		# if k > 1, update the old regression coefficients
		if( (k > 1) && (any(oldInNewLogical)) ) {
			betaHat[oldActiveIndex, k] = betaHat[oldActiveIndex, k-1] + oldBetaHatIncrement
		}
		
		# update the new regression coefficient, first find the correct position
		if( length( match(oldCurrentList, currentList) ) > 0 ){
			addVariable = currentList[-match(oldCurrentList, currentList)]
			} else {
				addVariable = currentList
			}
		
		addIndex = match(addVariable, attributes(myXStar)$names)
		betaHat[addIndex, k] = newBetaHatIncrement

		# lasso modification step: 

		if(k > 1){
			if(lassoFlag == TRUE){
				removeIndexTotal = match(cStar, newcVec)
				removeIndex = match(attributes(myXStar)$names[removeIndexTotal], currentList)
				
				potentialList = c(potentialList, currentList[removeIndex])
				currentList = currentList[-removeIndex]
				activeIndex = match(currentList, attributes(myXStar)$names)
				inactiveIndex = match(potentialList, attributes(myXStar)$names)
				

			}

		}

		k = k + 1
		np = np + 1
		if(lassoFlag){
			np = np - 1
		}

		if((np == min(p, n-1)-1) && (lassoFlag == TRUE) ){
			break; 
		}

	}

	# the np = min(p, n-1) step

	rp = myY - mu[, dim(mu)[2]]
	if( (length(potentialList) > 0) &&  (np < min(p, n-1)) ){
		XPotential = myXStar[potentialList]
		j = which.max( abs( cor( rp, XPotential ) ) )

		oldCurrentList = currentList
		oldPotentialList = potentialList
		currentList = c(currentList, potentialList[j])
		potentialList = potentialList[-j]  
		} 




	AkCorSign = as.vector(sign( t(as.matrix(myXStar[currentList])) %*% rp ))

	designMat = as.matrix(myXStar[currentList]) %*% diag(AkCorSign)

	# if introducing this parameter makes the design matrix singular, roll back to previous 
	if(det( t(designMat) %*% designMat ) == 0){
		currentList = oldCurrentList
		potentialList = oldPotentialList
	}
	AkCorSign = as.vector(sign( t(as.matrix(myXStar[currentList])) %*% rp ))
	designMat = as.matrix(myXStar[currentList]) %*% diag(AkCorSign)

	activeIndex = match(currentList, attributes(myXStar)$names)
	inactiveIndex = match(potentialList, attributes(myXStar)$names)
	oldActiveIndex = match(oldCurrentList, attributes(myXStar)$names)

	updateBeta = solve( (t(designMat) %*% designMat), (t(designMat) %*% (myY - mean(myY))) )
	updateMu = designMat %*% updateBeta + mean(myY)

	# this step is to restore the correct sign of betaHat
	updateBeta = AkCorSign * updateBeta

	if(dim(betaHat)[2] >= min(p, n-1)){
		betaHat = cbind(betaHat, rep(0, p))
	}
	betaHat[activeIndex, dim(betaHat)[2]] = updateBeta
	mu = cbind(mu, updateMu)


	# do some formatting for the output
	betaHat = cbind(rep(0, p), betaHat)
	betaHat = t(betaHat)
	

	# the rss
	rss = apply(mu, 2, function(xx){
		return( t(myY - xx) %*% (myY - xx) )
		})


	# create a lassoClass object to store the result for returning 
	lassoResult = new("lassoClass", betaHatClass = betaHat, muClass = mu, rssClass = rss)

	return(lassoResult)

}