\name{plot.StabilityReport}
\alias{plot.StabilityReport}
\title{Display the stability of a clustering as a heat map plot.}
\usage{
  \method{plot}{StabilityReport} (x, classes = NULL,
    class_colors = NULL, sort.clusters = 0, ...)
}
\arguments{
  \item{x}{A StabilityReport object, as given by an output
  of perturbationStability.}

  \item{classes}{Auxiliary class labels for the data
  points, possibly from known classes or other clusterings.
  The classes must be integers in 1,...,L.  If NULL, this
  column is not plotted.}

  \item{class_colors}{Colors to use when plotting the
  auxiliary class labels. If the given classes are in
  1,...,L, it must be a list of at least L colors. If NULL,
  \code{RColorBrewer} is used to choose representative
  colors. Ignored if \code{classes} is \code{NULL}.}

  \item{sort.clusters}{Whether to sort the clusters in the
  stability map image for aesthetic reasons.  0 (default)
  means to not reorder them, 1 orders them by cluster size,
  and 2 orders them by average stability.}

  \item{...}{optional arguments passed to internal
  functions}
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
\examples{
## Generate a fake dataset with 3 clusters
cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
cl.size <- 100
X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
            rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
            rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))
dists  <- t(apply(X, 1, function(mu) {sqrt(rowSums((cen - mu)^2))}))
labels <- c(rep(1,100), rep(2,100), rep(3,100))

## Apply to just the distance matrix
stability1 <- perturbationStability(dists)

## Ways to display information
print(stability1)
summary(stability1)
plot(stability1, classes=labels)

## Add in our labels
cl <- list(dists = dists, labels = labels)
stability2 <- perturbationStability(cl)

print(stability2)
summary(stability2)
plot(stability2, classes=labels)
}
\seealso{
  \code{\link{easystab}}
}

