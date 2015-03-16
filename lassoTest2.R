# lassoTest 2

source("lasso.R")

library(lars)
data(diabetes)
myX = diabetes$x
myY = diabetes$y

myX = apply(myX, 2, function(xx){
	return( ((xx - mean(xx))/sd(xx)))
	})


mylasso1 = lars(myX, myY, "lasso")

mylar1 = lars(myX, myY, "lar")

mylasso2 = lasso(myY, myX)

dev.new()
plot(mylasso1, breaks = FALSE)
#dev.off()

dev.new()
#pdf("mylasso.pdf")
plotCoef(mylasso2)
#dev.off()

print(coef(mylasso1))

print(betaHatClass(mylasso2))