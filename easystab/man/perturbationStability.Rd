\name{perturbationStability}
\alias{perturbationStability}
\title{calculate perturbation stability of a list of clusterings}
\description{calculate perturbation stability of a list of clusterings}
\usage{
perturbationStability(clusterings, n_baselines = 32, seed = 0, use_permutations = FALSE, by_dimension = FALSE, Kmap_mode = 2, theta = NULL)
}
\arguments{
  \item{clusterings}{A list of clustering. Each clustering is a list
    with a distance to clustering centroid matrix}

  \item{n_baselines}{n_baselines}

  \item{seed}{seed}

  \item{use_permutations}{use_permutations}

  \item{by_dimension}{by_dimension}

  \item{Kmap_mode}{Kmap_mode}

  \item{theta}{If null, an optimized theta will be calcuated within function}
}
\value{
  A list of new cluserings. Each clustering has extra components
  including \code{sorted_stability_matrix},
  \code{sorted_stability_matrix_index_map},
  \code{sorted_stability_matrix_cluster_map}, \code{opt_theta}
}