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

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

dx <- dist(scale(X))
hc <- hclust(dx)

cl_list <- from.hclust(dx, hc)
stability_sequence <- perturbationStability(cl_list)

# Information about the stability sequence
print(stability_sequence)
summary(stability_sequence)

# Plot the stability sequence
plot(stability_sequence)

# Plot the most stable clustering.
plot(stability_sequence$best, classes = iris[,"Species"])
}
\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}},
  \code{\link{clusterings_from_kmeans}}
}

