\name{plotStabilityMap}

\alias{plotStabilityMap}

\title{Represent the stability of a clustering as a heat map plot.}

\description{Plots the stability of a clustering as a heat map plot,
  showing the relative stability of the different clusters and the
  overall behavior of the clustering.  The input is taken as a single
  clustering analysis as given by perturbationStability.
}

\usage{
plotStabilityMap(clustering)
}
\arguments{
  \item{clustering}{A clustering with stability analysis results
    attached, as given by an output of perturbationStability.}
}

\seealso{
  \code{\link{easystab}},   \code{\link{plotLabeledStabilityMap}}
} 