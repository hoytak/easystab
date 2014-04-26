\name{from.hclust}
\alias{from.hclust}
\title{Adapts the output of \code{hclust} for input into
\code{perturbationStability}.}
\usage{
  from.hclust(dx, hc, k = 1:10, method = "average")
}
\arguments{
  \item{dx}{Distance matrix as produced by \code{dists},
  giving the point-to-point distances.}

  \item{hc}{Hierarchical clustering as produced by
  \code{hclust}.}

  \item{k}{A list giving the numbers of clusters to cut the
  tree at; this is passed to \code{\link{cutree}}.
  Defaults to 1:10.}

  \item{method}{Method used to calculate the
  point-to-cluster distances from the point-to-point
  distance matrix \code{dx} given.  Currently, the two
  supported methods are "average", which takes the average
  of the distances between the given point and all the
  points in the cluster (similar to average linkage), and
  "median", which uses the median distance to the points in
  the cluster.}
}
\value{
  A list of clusterings suitable for use with
  \code{perturbationStability}.
}
\description{
  Adapts the output of \code{hclust} for use with
  \code{perturbationStability} to give more information
  about the behavior of the hierarchical clustering tree.
}
\examples{
############################################################
## Interfacing with the hierarchical clustering method
library(easystab)

## Generate a fake dataset with 3 clusters
cen <- matrix(c(0,-2,1,3,-3,1), ncol=2, byrow=TRUE)
cl.size <- 100
X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
           rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
           rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))

dx <- dist(X)
hc <- hclust(dx)
cl_list <- from.hclust(dx,hc)

stability_collection <- perturbationStability(cl_list)

## Information about the stability sequence
print(stability_collection)
summary(stability_collection)

## Plot the stability sequence
plot(stability_collection)

############################################################
## A more detailed example using the UCI Wisconsin breast cancer dataset.
library(mlbench)

# Load and cluster the Breast Cancer dataset using correlation distance.
data(BreastCancer)

bcdata <- na.omit(BreastCancer)

## Use 1 - (x \%*\% y) / (|x|_2 |y|_2) to compute divergence
X <- data.matrix(bcdata[,-c(1,11)])
Y <- X \%*\% t(X)
Ynorm <- diag(diag(Y)^(-1/2))
dx <- as.dist(1 - Ynorm \%*\% Y \%*\% Ynorm)
hc <- hclust(dx, method="complete")

cl_list <- from.hclust(dx, hc, method = "median")
stability_collection <- perturbationStability(cl_list)

# Information about the stability sequence
print(stability_collection)
summary(stability_collection)

layout(matrix(1:2, nrow=1, ncol=2))
plot(stability_collection)
plot(stability_collection$best, classes = bcdata[,11])
}
\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}},
  \code{\link{from.kmeans}}
}

