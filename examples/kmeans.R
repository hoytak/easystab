library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

km_list = list()

for(k in 1:10) {
  km_list[[k]] = kmeans(X, k)
}

stability_sequence <- perturbationStability(kmeans_stability(X, km_list))

plot(stability_sequence)
# plot(stability_sequence[[5]], with_label = TRUE, classes = iris[,"Species"])
