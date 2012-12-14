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
library(easystab)
data(sample)
pl<-perturbationStability(clusterings)
res<-estimateK(pl)
l<-pl[[res$estimate_index]]
plotStabilitySequence(clusterings)
}
