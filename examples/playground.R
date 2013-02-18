## Toy example using a simple mixture model
library(easystab)

cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)

cl.size <- 100
## X is 300 x 2 matrix of 2d points, 100 from each of 3 components
X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]),
                   rnorm(cl.size,mean=cen[[1,2]])),
            rbind(rnorm(cl.size,mean=cen[[2,1]]),
                  rnorm(cl.size,mean=cen[[2,2]])),
            rbind(rnorm(cl.size,mean=cen[[3,1]]),
                  rnorm(cl.size,mean=cen[[3,2]]))))

dists  <- t(apply(X, 1, function(mu) {sqrt(rowSums((cen - mu)^2))}))
labels <- c(rep(1,cl.size), rep(2,cl.size), rep(3,cl.size))

## Apply to just the distance matrix
stability1 <- perturbationStability(dists)

## Ways to display information
print(stability1)
summary(stability1)
plot(stability1, classes=labels)

## Add in our labels
cl <- list(dists = dists, labels = labels)
stability2 <- perturbationStability(cl)

print(stability2)
summary(stability2)
plot(stability2, classes=labels)

## Now try several numbers of clusters using kmeans
km_list = list()

for(k in 1:8)
 km_list[[k]] = kmeans(X, k, iter.max=50, nstart=100)

cl_list <- from.kmeans(X, km_list)
stability_sequence <- perturbationStability(cl_list)

print(stability_sequence)
summary(stability_sequence)
plot(stability_sequence)

layout(matrix(1:9, ncol = 3, byrow=TRUE))

for(i in 1:9) {
  t <- (i - 1)/4
  Z <- make2dStabilityImage(cen, theta=t, buffer=2)
  image(Z$x, Z$y, Z$stability, main = sprintf("Theta = %1.2f.", t),
        xlab = "x", ylab="y")
}


## Toy example using a simple mixture model
cen <- matrix(c(0,0,1,2,-1,0.5), ncol=2)

## X is 300 x 2 matrix of 2d points, 100 from each of 3 components
X <- t(cbind(rbind(rnorm(100,mean=cen[[1,1]],sd=1),rnorm(100,mean=cen[[1,2]],sd=1)),
            rbind(rnorm(100,mean=cen[[2,1]],sd=2),rnorm(100,mean=cen[[2,2]],sd=2)),
            rbind(rnorm(100,mean=cen[[3,1]],sd=1),rnorm(100,mean=cen[[3,2]],sd=1))))

dists  <- t(apply(X, 1, function(mu) {sqrt(rowSums((cen - mu)^2))}))
labels <- c(rep(1,100), rep(2,100), rep(3,100))

## Apply to just the distance matrix
stability1 <- perturbationStability(dists)

## Ways to display information
print(stability1)
summary(stability1)
plot(stability1, classes=labels)

## Add in our labels
cl <- list(dists = dists, labels = labels)
stability2 <- perturbationStability(cl)

print(stability2)
summary(stability2)
plot(stability2, classes=labels)

## Now try several numbers of clusters using kmeans
km_list <- list()

for(k in 1:8)
 km_list[[k]] <- kmeans(X, k, iter.max=50, nstart=100)

cl_list <- from.kmeans(X, km_list)

stability_sequence <- perturbationStability(cl_list)
plot(stability_sequence)

## Now plot each K with multiple runs of the clustering function.
## Now try several numbers of clusters using kmeans
km_list <- list()

km_list <- lapply(0:23, function(k) { kmeans(X, 1 + (k %% 8))})
stability_sequence <- perturbationStability(from.kmeans(X, km_list))
plot(stability_sequence)



## Now try a range of numbers of clusters using kmeans
km_list1 <- lapply(1:6, function(k) { kmeans(X, k, iter.max=50, nstart=100)})
stabilities1 <- perturbationStability(from.kmeans(X, km_list1))

plot(stabilities1)

## Now plot each K with multiple runs of the clustering function.
## Now try several numbers of clusters using kmeans
km_list2 <- lapply(0:17, function(k) { kmeans(X, 1 + (k %% 6))})
stabilities2 <- perturbationStability(from.kmeans(X, km_list2))

plot(stabilities2)

## Plot the same thing, except without grouping by number of clusters
plot(stabilities2, sort=FALSE)

## If two clusterings have the same number of clusters, plot only the
## most stable one.
plot(stabilities2, prune=TRUE, sort=FALSE)

## Name the best one
stabilities2[[stabilities2$best.index]]$name <- "BEST!!!"
plot(stabilities2)
