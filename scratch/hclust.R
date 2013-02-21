library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

dx <- dist(scale(X))
hc <- hclust(dx)

cl_list <- hclust_stability(dx, hc)

stability_sequence <- perturbationStability(cl_list)

print(stability_sequence)

summary(stability_sequence)

plot(stability_sequence)

# plot(stability_sequence[[3]], with_label = TRUE, classes = iris[,"Species"])
