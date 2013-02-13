\name{plot.StabilityReport}
\alias{plot.StabilityReport}
\title{Display the stability of a clustering as a heat map plot.}
\usage{
  plot.StabilityReport(clustering, classes = NULL,
    class_colors = NULL)
}
\arguments{
  \item{clustering}{A StabilityReport object, as given by
  an output of perturbationStability.}

  \item{classes}{Auxiliary class labels for the data
  points, possibly from known classes or other clusterings.
  The classes must be integers in 1,...,L.  If NULL, this
  column is not plotted.}

  \item{class_colors}{Colors to use when plotting the
  auxiliary class labels. If the given classes are in
  1,...,L, it must be a list of at least L colors. If NULL,
  \code{RColorBrewer} is used to choose representative
  colors. Ignored if \code{classes} is \code{NULL}.}
}
\description{
  Plots the stability of a clustering as a heat map plot,
  showing the relative stability of the different clusters,
  the data points, and the overall behavior of the
  clustering.  The input is taken as a single clustering
  analysis as given by perturbationStability.
}
\details{
  If \code{classes} are supplied (possibly from known
  classes or from another clustering) version, they are
  plotted alongside the heatmap plot, with class membership
  indexed by color.
}
\seealso{
  \code{\link{easystab}}
}

