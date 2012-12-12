\name{perturbationStability}
\alias{perturbationStability}
\title{Calculate perturbation stability of a list of clusterings}
\description{Calculate the perturbation stability of a list of
  clusterings. 

}
\usage{
perturbationStability(clusterings, n_baselines = 32, seed = 0, use_permutations = FALSE, by_dimension = FALSE, Kmap_mode = 0, theta = NULL)
}
\arguments{
  \item{clusterings}{A list of clustering. Each clustering contains a
    distance to centroid or distance to clustering matrix.}

  \item{n_baselines}{n_baselines}  The number of random baseline
  matrices to use in computing the distance 

  \item{seed}{seed}

  \item{use_permutations}{use_permutations}

  \item{by_dimension}{by_dimension}

  \item{Kmap_mode}{Kmap_mode} Determines how the clusters in the cluster
   stability plot are organized. Kmap_mode = 0 (default) does not
   reorder the clusters from the original label map. Kmap_mode = 1 sorts
   them by cluster size, and Kmap_mode = 2 sorts them by cluster 
   are sorted 

  \item{theta}{If null, an optimized theta will be calcuated within function}
}
\value{
  A list of new cluserings. Each clustering has extra components
  including \code{sorted_stability_matrix},
  \code{sorted_stability_matrix_index_map},
  \code{sorted_stability_matrix_cluster_map}, \code{opt_theta}
}