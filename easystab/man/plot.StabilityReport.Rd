\name{plot.StabilityReport}

\alias{plot.StabilityReport}

\title{Plots the stability report of a clustering as a heat map plot.}

\description{Plots the stability of a clustering as a heat map plot,
  showing the relative stability of the different clusters and the
  overall behavior of the clustering.  The input is taken as a single
  clustering analysis as given by perturbationStability.

  If \code{classes} are supplied (possibly from known classes or from another
  clustering) version, they are plotted alongside the heatmap plot, with
  class membership indexed by color.
}

\usage{
plot.StabilityReport(clustering, with_label = FALSE, classes = NULL, class_colors = NULL)
}
\arguments{
  \item{clustering}{A clustering with stability analysis results
    attached, as given by an output of perturbationStability.}
  
  \item{with_label}{option to enable cluster label.}  

  \item{classes}{Auxiliary class labels for the data points, possibly from
    known classes or other clusterings. The classes must integers in
    1,...,L.  If NULL, clustering's own labels will be used. Ignored
    when with_label is FALSE.}  

  \item{class_colors}{Colors to use when plotting the auxiliary class
    labels.  If the given classes are in 1,...,L, it must be a list of
    at least L colors.  If NULL, \code{RColorBrewer} is used to choose
    representative colors.  }
}

\seealso{
  \code{\link{easystab}}
} 