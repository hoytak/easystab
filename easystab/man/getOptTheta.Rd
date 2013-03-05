\name{getOptTheta}
\alias{getOptTheta}
\title{Calculates the optimal prior parameter}
\usage{
  getOptTheta(clusterings, seed = 0, n_baselines = 25)
}
\arguments{
  \item{clusterings}{A single clustering or a list of
  clusterings. Each clustering of \code{n} data points into
  \code{K} clusters is specified primarily by matrix giving
  point to cluster distances.  Specifically, clustering
  must contain an of \code{n} by \code{K} distance matrix
  giving the point to cluster distance (\code{K} can be
  different across clusterings). Optionally, an array of
  \code{n} integer labels \code{labels} in \code{1,...,K}
  is expected; if not present, a warning is given and the
  labels are computed according to the minimum
  point-to-cluster distance.}

  \item{seed}{Random seed used for generating the baseline
  stability matrices.}

  \item{n_baselines}{The number of random baseline matrices
  to use in computing the stability scores.  Increase this
  number to get more accuracy at the expense of speed.}
}
\description{
  Calculate the optimal prior parameter theta by maximizing
  the difference in overall stability aganst the baseline
  distributions. The theta parameter indexes the strength
  of the perturbations, with smaller values translating
  into stronger perturbations.
}
\examples{
##################################################
## These examples produce exactly the same results as those in
## perturbationStability.

library(easystab)
## Generate a fake dataset with 3 clusters

cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
cl.size <- 100
X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]),
                  rnorm(cl.size,mean=cen[[1,2]])),
                        rbind(rnorm(cl.size,mean=cen[[2,1]]),
                              rnorm(cl.size,mean=cen[[2,2]])),
                        rbind(rnorm(cl.size,mean=cen[[3,1]]),
                              rnorm(cl.size,mean=cen[[3,2]]))))
dists  <- t(apply(X, 1, function(mu) {sqrt(rowSums((cen - mu)^2))}))
labels <- c(rep(1,100), rep(2,100), rep(3,100))

## Takes same input as
theta <- getOptTheta(dists)

## Apply to just the distance matrix
stability1 <- perturbationStability(dists, theta = theta)

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
km_list <- lapply(1:8, function(k) { kmeans(X, k, iter.max=20, nstart=30)})
cl_list <- from.kmeans(X, km_list)
stability_collection <- perturbationStability(cl_list)

print(stability_collection)
summary(stability_collection)
plot(stability_collection)
}
\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}}
}

