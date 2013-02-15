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
#'@param clustering StabilityReport object. A clustering with stability
#'analysis results attached, as given by an output of perturbationStability.
#'@param with_label option to enable cluster label.
#'@param classes Auxiliary class labels for the data points, possibly from
#'known classes or other clusterings. The classes must integers in 1,...,L.  If
#'NULL, clustering's own labels will be used. Ignored when with_label is FALSE.
#'@param class_colors Colors to use when plotting the auxiliary class labels.
#'If the given classes are in 1,...,L, it must be a list of at least L colors.
#'If NULL, \code{RColorBrewer} is used to choose representative colors.
#'@param clusterings StabilitySequence object. The output of
#'\code{perturbationStablity} -- a list of clusters with perturbation stability
#'analyses.
#'@seealso \code{\link{perturbationStability}},
#'\code{\link{plotStabilityImage}}
#'@examples
#'
#'
#'## example with built in sample data
#'## clusterings in sample data already have 'dist' matrix included
#'
#'library(easystab)
#'
#'data(sample)
#'pl<-perturbationStability(clusterings)
#'res<-estimateK(pl)
#'
#'l<-pl[[res$estimated_index]]
#'
#'plot(pl)
#'plot(l)
#'
#'## example with hclust function on iris data set
#'
#'X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]
#'
#'dx <- dist(scale(X))
#'hc <- hclust(dx)
#'cl_list <- hclust_stability(dx, hc)
#'
#'stability_sequence <- perturbationStability(cl_list)
#'
#'print(stability_sequence)
#'
#'summary(stability_sequence)
#'
#'## example with kmeans function on iris data set
#'
#'X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]
#'
#'km_list = list()
#'
#'for(k in 1:10) {
#'  km_list[[k]] = kmeans(X, k)
#'}
#'
#'stability_sequence <- perturbationStability(kmeans_stability(X, km_list))
#'
#'# plots the sequence
#'plot(stability_sequence)
#'
#'# plots the stability map
#'plot(stability_sequence[[3]], with_label = TRUE, classes = iris[,"Species"])
#'
#'
NULL
