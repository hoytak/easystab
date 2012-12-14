\name{plotLabeledStabilityMap}

\alias{plotLabeledStabilityMap}

\title{Plot the stability heat map, with labels attached, based on clustering stability
  analysis results.}

\description{Plots the stability of a clustering as a heat map plot,
  showing the relative stability of the different clusters and the
  overall behavior of the clustering.  The input is taken as a single
  clustering analysis as given by perturbationStability.  In this
  version, supplied class labels from another source (possibly from
  known classes or from another clustering), are plotted alongside

}
\usage{
plotLabeledStabilityMap(clustering, labels, label_colors = NULL)
}
\arguments{
  \item{clustering}{A clustering with stability analysis results
    attached.}
  \item{labels}{Auxiliary labels for the data points, possibly from
    known classes or other clusterings. The labels must integers in
    1,...,L.}  

  \item{label_colors}{Colors to use when plotting the auxiliary labels.
    If the labels are in 1,...,L, it must be a list of L colors.   If
    NULL, \code{RColorBrewer} is used to choose representative colors.
  }

\value{
  \item{colors}{The colors used to plot the supplied labels.}
}

\seealso{
  \code{\link{easystab}}, \code{\link{perturbationStability}},
  \code{\link{plotStabilityMap}}
}

