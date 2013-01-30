\name{easystab-package}

\alias{easystab-package}

\alias{easystab}

\docType{package}

\title{ Clustering stability analysis using Bayesian perturbations.  }

\description{ A collection of functions for analyzing the stability of
  one or more clusterings.  The technique used is given in [[PAPER]].
  Functions are given to: assess the behavior of a clustering under
  perturbations, estimate the correct number of clusters in a dataset,
  assess the relative stability of individual clusters within a
  dataset, and interface with \code{hclust} and \code{kmeans}. }


\details{
  \tabular{ll}{
    Package: \tab easystab\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2012-11-28\cr
    License: \tab GPL (>= 2)\cr
    LazyLoad: \tab yes\cr
  }
}

\examples{

## example with built in sample data
## clusterings in sample data already have 'dist' matrix included

library(easystab)

data(sample)
pl<-perturbationStability(clusterings)
res<-estimateK(pl)

l<-pl[[res$estimated_index]]

plotStabilitySequence(pl)
plotStabililityMap(l)

## example with hclust function on iris data set

library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]


dx <- dist(scale(X))
hc <- hclust(dx)
cl_list <- hclust_stability(dx, hc)

stab_list <- perturbationStability(cl_list)

plotStabilitySequence(stab_list)

## example with kmeans function on iris data set

library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

km_list = list()

for(k in 1:10) {
  km_list[[k]] = kmeans(X, k)
}

stab_list <- perturbationStability(kmeans_stability(X, km_list))

plotStabilitySequence(stab_list)

plotStabilityMap(stab_list[[3]], with_label = TRUE, classes = iris[,"Species"])

}
