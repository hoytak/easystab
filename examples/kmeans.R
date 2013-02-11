library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

km_list = list()

for(k in 1:8) {
  km_list[[k]] = kmeans(X, k, iter.max=50, nstart=100)
}

stability_sequence <- perturbationStability(kmeans_stability(X, km_list))

plot(stability_sequence$best, classes = iris[,"Species"])
# plot(stability_sequence[[5]], with_label = TRUE, classes = iris[,"Species"])
