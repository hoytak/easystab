\name{hclust_stability}

\alias{hclust_stability}

\title{Adapts the output of \code{hclust} for input into \code{perturbationStability}.}

\description{Adapts the output of \code{hclust} for use with
  \code{perturbationStability} to give more information about the
  behavior of the hierarchical clustering tree. }

\usage{
hclust_stability(dx, hc, clsnum_min = 1, clsnum_max = 10, method = "average")
}

\arguments{
  \item{dx}{Distance matrix as produced by \code{dists}, giving
    the point-to-point distances.}

  \item{hc}{Hierarchical clustering as produced by \code{hclust}.}
  
  \item{clsnum_min}{Minimum cluster number (default 1).}
  
  \item{clsnum_max}{Maximum cluster number (default 10).}
  
  \item{method}{ Method used to calculate the point-to-cluster distances
    from the point-to-point distance matrix \code{dx} given.  Currently,
    the two supported methods are "average", which takes the average of
    the distances between the given point and all the points in the
    cluster (similar to average linkage), and "median", which uses
    median distance to the points in the cluster. }
  
}

\value{A list of clusterings suitable for use with \code{perturbationStability}.}

\seealso{
  \code{\link{easystab}}, \code{\link{perturbationStability}},
  \code{\link{clusterings_from_kmeans}}
}

