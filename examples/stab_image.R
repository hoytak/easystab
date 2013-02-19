library(easystab)

cen <- matrix(c(0,-1,-2,2,0,1, 3,3), nrow=4, byrow=TRUE)

Z <- make2dStabilityImage(cen, theta=2, size=c(502,502), buffer=2)

image(Z$x, Z$y, Z$stability)
