\name{perturbationStability}

\alias{perturbationStability}

\title{Calculate perturbation stability of a list of clusterings}

\description{Calculate the perturbation stability of a list of
  clusterings. }

\usage{ perturbationStability(clusterings, n_baselines = 32, seed = 0,
Kmap_mode = 0, theta = NULL) }

\arguments{ 

  \item{clusterings}{A list of clusterings. For each clustering of
\code{n} points into \code{K} clusters (\code{K} can be different across
clusterings), an \code{n} by \code{K} distance to centroid or distance
to clustering matrix \code{dists} and array of \code{n} integer labels
\code{labels} in \code{1,...,K} is expected.}

  \item{n_baselines}{n_baselines}  The number of random baseline
  matrices to use in computing the stability scores.  Increase this
  number to get more accuracy 

  \item{seed}{seed} The random seed used for sampling the baseline
  matrices.

  \item{Kmap_mode}{Kmap_mode} Determines how the clusters in the cluster
   stability heatmap are organized. Kmap_mode = 0 (default) does not
   reorder the clusters from the original label map. Kmap_mode = 1 sorts
   them by cluster size, and Kmap_mode = 2 sorts them by individual
   cluster stability.

   \item{theta}{The non-negative rate parameter passed to the shifted exponential prior on
     the perturbations.  The parameter indexes the strength of the
     perturbations, with smaller values translating into stronger
     perturbations.  If NULL, theta is chosen by optimizing the overall
     stability against the baseline distributions.}
 }
 
 \value{
   Returns a list of clusterings and stability information matching the
  supplied list of clusterings.  Each object in this list contains the
  original \code{dists} and \code{labels}, the scalar stability score
  \code{stability_score}, the empirical collection of stability scores
  \code{scores}, the theta parameter used \code{theta}, and the components of
  the sorted stability matrix used for plotting the behavior of the
  clustering.
}

