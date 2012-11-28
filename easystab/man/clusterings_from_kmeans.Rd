\name{clusterings_from_kmeans}
\alias{clusterings_from_kmeans}
\title{perturbationStability function clustering input adaptor for kmeans}
\description{Run kmeans clustering analysis and return list of
  clusterings as perturbationStability function input}
\usage{
clusterings_from_kmeans(x, clsnum_min = 2, clsnum_max=10)
}
\arguments{
  \item{x}{matrix or data frame object}
  \item{clsnum_min}{minimum cluster number. default to 2}
  \item{clsnum_max}{maximum cluster number. default to 10}
}
\seealso{
  \code{\link{easystab}}, \code{\link{perturbationStability}},
  \code{\link{clusterings_from_hclust}}
}

