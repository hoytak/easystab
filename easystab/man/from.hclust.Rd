\name{from.hclust}
\alias{from.hclust}
\title{Adapts the output of \code{hclust} for input into
\code{perturbationStability}.}
\usage{
  from.hclust(dx, hc, clsnum_min = 1, clsnum_max = 10,
    method = "average")
}
\arguments{
  \item{dx}{Distance matrix as produced by \code{dists},
  giving the point-to-point distances.}

  \item{hc}{Hierarchical clustering as produced by
  \code{hclust}.}

  \item{clsnum_min}{Minimum cluster number (default 1).}

  \item{clsnum_max}{Maximum cluster number (default 10).}

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
cen <- matrix(c(0,-20,10,20,-20,10), ncol=2, byrow=TRUE)
cl.size <- 100
X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
            rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
            rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))

dx <- dist(X)
hc <- hclust(dx)
cl_list <- from.hclust(dx,hc)

stability_collection <- perturbationStability(cl_list)

# Information about the stability sequence
print(stability_collection)
summary(stability_collection)

# Plot the stability sequence
plot(stability_collection)

## Ways to display information
print(stability1)
summary(stability1)
plot(stability1, classes=labels)










library(lsa)

breast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",sep=",")

X <- breast[,-c(1,2)]

dx <- as.dist(cosine(t(X)))
hc <- hclust(dx, "ward")

cl_list <- from.hclust(dist(scale(X)), hc)
stability_collection <- perturbationStability(cl_list)

# Information about the stability sequence
print(stability_collection)
summary(stability_collection)

# Plot the stability sequence
plot(stability_collection)

# Plot a map of the most stable clustering.
plot(stability_collection$best, classes = breast[,2])
}
\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}},
  \code{\link{}}
}

