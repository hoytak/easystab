library(easystab)
library(lsa)

breast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",sep=",")

X <- breast[,-c(1,2)]

dx <- as.dist(1 - cosine(t(X)))
hc <- hclust(dx)

cl_list <- from.hclust(dx, hc)
stability_collection <- perturbationStability(cl_list)

# Information about the stability sequence
layout(matrix(1:2, nrow=1, ncol=2))
print(stability_collection)
summary(stability_collection)

plot(stability_collection)

plot(stability_collection$best, classes = breast[,2])
