library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

km_list = list()

for(k in 1:10) {
  km_list[[k]] = kmeans(X, k)
}

stab_list <- perturbationStability(kmeans_stability(X, km_list))

plotStabilitySequence(stab_list)

plotStabilityMap(stab_list[[3]], with_label = TRUE, classes = iris[,"Species"])
