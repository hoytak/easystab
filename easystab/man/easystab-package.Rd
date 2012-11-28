\name{easystab-package}
\alias{easystab-package}
\alias{easystab}
\docType{package}
\title{
  Clustering perturbabtion stability analysis
}
\description{
  Package of functions for clustering perturbation stability analysis
}
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
