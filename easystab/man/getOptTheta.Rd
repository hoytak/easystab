\name{getOptTheta}
\alias{getOptTheta}
\title{Calculates the optimal prior parameter 

  \description{Calculate the optimal prior parameter theta by maximizing
    the difference in overall stability aganst the baseline
    distributions. The theta parameter indexes the strength of the
    perturbations, with smaller values translating into stronger
    perturbations. 
  }
\usage{
getOptTheta(clusterings, seed = 0, n_baselines = 32)
}
\arguments{

  \item{clusterings}{A list of clusters. A list of clusterings. For each
clustering of \code{n} points into \code{K} clusters (\code{K} can be
different across clusterings), an \code{n} by \code{K} distance to
centroid or distance to clustering matrix \code{dists} and array of
\code{n} integer labels \code{labels} in \code{1,...,K} is expected.}

  \item{seed}{Random seed used for generating the baseline stability matrices.}

  \item{n_baselines}{ The number of random baseline matrices to use in
  computing the stability scores.  Increase this number to get more
  accuracy in distinguishing between clusterings. }
}

\seealso{
  \code{\link{easystab}}, \code{\link{perturbationStability}}
}

