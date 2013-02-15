\name{plot.StabilityCollection}
\alias{plot.StabilityCollection}
\title{Plot the stability scores produced by perturbationStability as a sequence of
box plots.}
\usage{
  \method{plot}{StabilityCollection} (clusterings,
    sort = TRUE, prune = FALSE, label.indices = NULL,
    color.best = "red", ...)
}
\arguments{
  \item{clusterings}{The output of
  \code{perturbationStablity} -- a list of clusters with
  perturbation stability analyses. Additionally:

  Set the \code{label} attribute of a specific clustering
  in order to change the corresponding label on the box
  plots. For example, \code{clusterings[[5]]$label <-
  "Clust5"} sets the displayed label of that clustering,
  overriding the generated labels.

  Set the \code{color} attribute of a specific clustering
  in order to change the color of boxplot.  The default is
  to color the "best" one red, and the rest black (see
  \code{color.best} below).  This overrides this behavior.}

  \item{sort}{Whether to sort the results in ascending
  order by the number of clusters in the data, then by
  stability scores within the clusters.}

  \item{prune}{If sort is TRUE, and multiple clusterings
  are given for a specific number of clusters, then show
  only the most stable one from each group.  For example,
  if there were three clusterings in the collection that
  had 5 clusters, only the most stable of those three would
  be displayed.}

  \item{label.indices}{If \code{label.indices} is TRUE,
  then the original indices from \code{clusterings} is
  included in the label for each box plot; if FALSE, they
  are not included.  If \code{label.indices} is NULL
  (default), then they are included only if items in the
  graph are reordered.  Note that setting the \code{label}
  attribute on the clusterings input overrides this.}

  \item{...}{Additional parameters passed to the boxplot
  function. See \code{\link{boxplot}} for more
  information.}

  \item{color.best}{Color the best clustering with this
  cluster. Ignored if NULL.}
}
\description{
  Summary display of the output of the
  perturbationStability function. Plot the stability scores
  produced by perturbationStability as a sequence of box
  plots.
}

