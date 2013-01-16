library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]


dx <- dist(scale(X))
hc <- hclust(dx)
cl_list <- hclust_stability(dx, hc)

stab_list <- perturbationStability(cl_list)

plotStabilitySequence(stab_list)
