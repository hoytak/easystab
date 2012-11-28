\name{estimateK}
\alias{estimateK}
\title{search cluster with maximum stability}
\description{parse result cluster list generated from
  pertubationStability function. find the cluster with maximum
  perturbation stability}
\usage{
estimateK(clusterings, p_value_threshold)
}
\arguments{
  \item{clusterings}{A list of clusters with perturbation stability
    analysis results. It is the result from perturbationStablity function}
  \item{p_value_threshold}{p_value_threshold. optional. default to 0.05}
}
\value{
  A list object with components estimate_K and estimate_index
}