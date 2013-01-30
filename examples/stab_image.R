library(easystab)
cen <- matrix(c(0.5,0,-1,0,0,1), nrow=3, byrow=TRUE)
plotStabilityImage(cen, 1.0, 1001, 1001, -2, 2, -2, 2)
