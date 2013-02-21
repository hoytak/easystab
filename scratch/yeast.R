library(easystab)

yeast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/yeast/yeast.data") 

X <- yeast[,-c(1,10)]
class_label <- as.numeric(yeast[,10])

# Sphere the data -- gives better results with 
X <- scale(X)
s <- svd(t(X))
X <- t(diag(1.0 / sqrt(s$d)) %*% t(s$u) %*% t(X))

km_list = list() 
for(k in 1:12) {
  km_list[[k]] = kmeans(X, k, iter.max = 50, nstart=50)
}

stability_sequence <- perturbationStability(kmeans_stability(X, km_list)) 

plot(stability_sequence)

