#'Clustering stability analysis using Bayesian perturbations.
#'
#'A collection of functions for analyzing the stability of one or more
#'clusterings.  The technique used is given in [[PAPER]].  Functions are given
#'to: assess the behavior of a clustering under perturbations, estimate the
#'correct number of clusters in a dataset, assess the relative stability of
#'individual clusters within a dataset, and interface with the \code{hclust}
#'and \code{kmeans} clustering packages.
#'
#'The exact method used is to perturb the cluster-to-point distances by scaling
#'them with a shifted exponential random variable, then computing the
#'probabilities of membership for each of the points under this perturbation.
#'This is compared against a set of randomly sampled bootstrapped baselines to
#'determine the final stability score.
#'
#'\tabular{ll}{ Package: \tab easystab\cr Type: \tab Package\cr Version: \tab
#'1.0\cr Date: \tab 2012-11-28\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
#'yes\cr }
#'
#'@name easystab-package
#'@aliases easystab-package easystab
#'@docType package
#'@seealso \code{\link{perturbationStability}}, \code{\link{make2dStabilityImage}},
#'\code{\link{plot.StabilityCollection}}, \code{\link{print.StabilityCollection}},
#'\code{\link{summary.StabilityCollection}}, \code{\link{plot.StabilityReport}},
#'\code{\link{print.StabilityReport}}, \code{\link{summary.StabilityReport}},
#'\code{\link{from.hclust}}, \code{\link{from.kmeans}}, \code{\link{getOptTheta}}
#'
#'@references Hoyt Koepke, Bertrand Clarke. to appear: Stastical Analysis and Data Mining.
#'
NULL
