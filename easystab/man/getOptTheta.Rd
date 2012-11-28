\name{getOptTheta}
\alias{getOptTheta}
\title{calculate optimized theta value for a list of clusters}
\description{calculate the minimum stability score of a list of
  clustering inputs}
\usage{
getOptTheta(clusterings, seed = 0, n_baselines = 32, use_permutation = FALSE , by_dimension = FALSE)
}
\arguments{
  \item{clusterings}{A list of clusters. Each cluster is expressed as
    a list with a component named as \code{dists}, which is the cluster point to centroid distance matrix}
  \item{seed}{seed, default to 0}
  \item{n_baselines}{n_baselines, default to 32}
  \item{use_permutation}{use_permutation switch, default to FALSE}
  \item{by_dimension}{by_dimension switch, default to FALSE}
}
\seealso{
  \code{\link{easystab}}, \code{\link{perturbationStability}}
}

