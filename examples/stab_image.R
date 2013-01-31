library(easystab)
cen <- matrix(c(1,0,-1,0,0,1), nrow=3, byrow=TRUE)
plotStabilityImage(cen, 0, 301, 301, -1.5, 1.5, -0.5, 2.5)
