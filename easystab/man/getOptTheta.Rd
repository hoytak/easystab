\name{getOptTheta}
\alias{getOptTheta}
\title{Calculates the optimal prior parameter}
\usage{
  getOptTheta(clusterings, seed = 0, n_baselines = 32)
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
\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}}
}

