# lassoRegressionClass Definition

source("interExtrapolate.R")

# define a lassoClass class to store the result of lasso
setClass("lassoClass", representation(betaHatClass = "matrix", muClass = "matrix", rssClass = "numeric"))

# accessor methods: 
#  betaHatClass():  gets the coefficients
#  muClass(): gets the fit 
#  rssClass(): gets the RSS


if(!isGeneric("betaHatClass")){
	if(is.function("betaHatClass")){
		fun = betaHatClass
	} else{
		fun = function(object) standardGeneric("betaHatClass")
		setGeneric("betaHatClass", fun)
	}
}

setMethod("betaHatClass", "lassoClass", function(object) object@betaHatClass)

if(!isGeneric("muClass")){
	if(is.function("muClass")){
		fun = muClass
	} else{
		fun = function(object) standardGeneric("muClass")
		setGeneric("muClass", fun)
	}
}

setMethod("muClass", "lassoClass", function(object) object@muClass)

if(!isGeneric("rssClass")){
	if(is.function("rssClass")){
		fun = rssClass
	} else{
		fun = function(object) standardGeneric("rssClass")
		setGeneric("rssClass", fun)
	}
}

setMethod("rssClass", "lassoClass", function(object) object@rssClass)

# additional class methods: 
#  predictClass: returns the fit for a newX predictor matrix at the step-th step of lasso
#  betaHatTClass: returns the coefficents when the complexity parameter t is newT (by interpolation)
#  predictTClass: returns the fit for a newX predictor matrix at when t = newT
#  tValueClass: returns the complexity parameters t at different steps of lasso

setGeneric("predictClass", function(object, newX, step) standardGeneric("predictClass"))

setMethod("predictClass", "lassoClass", function(object, newX, step = dim(object@betaHatClass)[1]) {

	predictYTemp = newX %*% object@betaHatClass[step, ] + object@muClass[1, 1]
	return(predictYTemp)

	})


setGeneric("betaHatTClass", function(object, newT) standardGeneric("betaHatTClass"))

setMethod("betaHatTClass", "lassoClass", function(object, newT=rowSums(abs(object@betaHatClass))[dim(object@betaHatClass)[1]] ){
	originalT = rowSums(abs(object@betaHatClass))
	p = dim(object@betaHatClass)[2]
	newBetaHat = rep(0, p)
	for(i in 1:p){
		originalBetaHatI = object@betaHatClass[, i]
		newBetaHat[i] = interExtrapolate(originalT, originalBetaHatI, newT)
	}
	return(newBetaHat)
	})




setGeneric("predictTClass", function(object, newX, newT) standardGeneric("predictTClass"))

setMethod("predictTClass", "lassoClass", function(object, newX, newT=rowSums(abs(object@betaHatClass))[dim(object@betaHatClass)[1]]) {

	originalT = rowSums(abs(object@betaHatClass))
	p = dim(object@betaHatClass)[2]
	newBetaHat = rep(0, p)
	for(i in 1:p){
		originalBetaHatI = object@betaHatClass[, i]
		newBetaHat[i] = interExtrapolate(originalT, originalBetaHatI, newT)
	}
	newX = as.matrix(newX)
	predictYTemp = newX %*% newBetaHat + object@muClass[1, 1]
	return(predictYTemp)

	})

setGeneric("tValuesClass", function(object) standardGeneric("tValuesClass"))

setMethod("tValuesClass", "lassoClass", function(object){
	originalT = rowSums(abs(object@betaHatClass))
	return(originalT)
	})


setGeneric("appErrorRateClass", function(object, myX, myY, newT = Inf) standardGeneric("appErrorRateClass"))

setMethod("appErrorRateClass", "lassoClass", function(object, myX, myY, newT = Inf){
	if(newT == Inf){
		newT = rowSums(abs(object@betaHatClass))[dim(object@betaHatClass)[1]]
	}

	myX = apply(myX, 2, function(xx){
		return( (xx - mean(xx, na.rm = TRUE))/sd(xx, na.rm = TRUE) )
		})

	originalT = rowSums(abs(object@betaHatClass))
	p = dim(object@betaHatClass)[2]
	newBetaHat = rep(0, p)
	for(i in 1:p){
		originalBetaHatI = object@betaHatClass[, i]
		newBetaHat[i] = interExtrapolate(originalT, originalBetaHatI, newT)
	}

	predictYTemp = myX %*% newBetaHat + object@muClass[1, 1]

	appErrorRate = sum( abs(myY - predictYTemp) )/length(myY)


	return(appErrorRate)

	})


setGeneric("noInfoErrorRateClass", function(object, myX, myY, newT = Inf) standardGeneric("noInfoErrorRateClass"))

setMethod("noInfoErrorRateClass", "lassoClass", function(object, myX, myY, newT = Inf){
	if(newT == Inf){
		newT = rowSums(abs(object@betaHatClass))[dim(object@betaHatClass)[1]]
	}

	myX = apply(myX, 2, function(xx){
		return( (xx - mean(xx, na.rm = TRUE))/sd(xx, na.rm = TRUE) )
		})

	originalT = rowSums(abs(object@betaHatClass))
	p = dim(object@betaHatClass)[2]
	newBetaHat = rep(0, p)
	for(i in 1:p){
		originalBetaHatI = object@betaHatClass[, i]
		newBetaHat[i] = interExtrapolate(originalT, originalBetaHatI, newT)
	}

	predictYTemp = myX %*% newBetaHat + object@muClass[1, 1]

	#appErrorRate = sum( (myY - predictYTemp)^2 )/length(myY)
	noInfoErrorRate = 0
	n = length(myY)
	for(tempI in 1:n){
		for(tempJ in 1:n){
			noInfoErrorRate = noInfoErrorRate +  abs( myY[tempI] - predictYTemp[tempJ] )
		}
	}
	noInfoErrorRate = noInfoErrorRate/(n^2)

	return(noInfoErrorRate)

	})

setGeneric("plotCoef", function(object) standardGeneric("plotCoef"))

setMethod("plotCoef", "lassoClass", function(object){
	step = dim(object@betaHatClass)[1]
	p = dim(object@betaHatClass)[2]
	tVal = tValuesClass(object)
	#dev.new()
	lineColors = rainbow(p)
	for(pIndex in 1:p){
		if(pIndex == 1){
			plot(tVal, object@betaHatClass[, pIndex], pch = ".", ylim= c(min(object@betaHatClass), max(object@betaHatClass) ))
			par(new = TRUE)
		}

		points(tVal, object@betaHatClass[, pIndex], pch = "*", col = lineColors[pIndex])
		lines(tVal, object@betaHatClass[, pIndex], col = lineColors[pIndex])
	}
	})


# a few additional methods that can be considered to implement: 
#  a method to compute apparent absolute predictive error rate 
#  a method to compute no information error rate
#  a method to plot those things


