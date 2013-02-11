############################################################
## Example using kmeans_stability on a single clustering

## Use yeast, X from previous example

## Works on a single clustering
km_cl <- kmeans(X, 8, iter.max = 50, nstart=50)
stability <- perturbationStability(kmeans_stability(X, km_cl))

## Plot the stability -- a single clustering, so displays it as a
## stability map plot.

plot(stability, classes = yeast[,10])
