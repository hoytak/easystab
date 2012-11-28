\name{clusterings_from_hclust}
\alias{clusterings_from_hclust}
\title{perturbationStability function clustering input adaptor for hclust}
\description{Run hclust clustering analysis and return list of
  clusterings as perturbationStability function input}
\usage{
clusterings_from_hclust(x, clsnum_min = 2, clsnum_max = 10, dist_method = "euclidean", method = "average")
}
\arguments{
  \item{x}{matrix or data frame object}
  \item{clsnum_min}{minimum cluster number. default to 2}
  \item{clsnum_max}{maximum cluster number. default to 10}
  \item{dist_method}{distance calculation method. default to
    \code{euclidean}}
  \item{method}{agglomeration method used in hclust. currently support
    \code{average} and \code{median} only}
}
\seealso{
  \code{\link{easystab}}, \code{\link{perturbationStability}},
  \code{\link{clusterings_from_kmeans}}
}

