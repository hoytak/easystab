\name{estimateK}
\alias{estimateK}

\title{Estimates the number of clusters.}
  
\description{

  Estimates the number of clusters in a dataset from a list of
clusterings and stability information as produced by
\code{perturbationStablity}.  The cluster is chosen as the one with the
fewest number of clusters that is statistically indistinguishable from
the most stable clustering with p-value \code{p_value}.

\usage{
estimateK(clusterings, p_value = 0.05)
}

\arguments{
  \item{clusterings}{The output of \code{perturbationStablity}
-- a list of clusters with perturbation stability analyses. }

  \item{p_value}{p-value requirement for saying that one clustering is
    more stable than another.}
}

\value{
  A list object with components estimated_K and estimated_index.
}

