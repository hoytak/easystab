\name{make2dStabilityImage}
\alias{make2dStabilityImage}
\title{Creates an image of the relative regions of stability for a 2d clustering.}
\usage{
  make2dStabilityImage(centroids, theta = 1, bounds = NULL,
    size = c(500, 500), buffer = 0.25)
}
\arguments{
  \item{centroids}{Array of 2D centroid points, given as a
  K by 2 array or matrix.}

  \item{theta}{The rate parameter passed to the shifted
  exponential prior on the perturbations.  \code{theta}
  must be non-negative; a warning is issued if \code{theta
  < 0}. The parameter indexes the strength of the
  perturbations, with smaller values translating into
  stronger perturbations.  If NULL, theta is chosen by
  optimizing the overall stability against the baseline
  distributions as in \code{\link{getOptTheta}}.}

  \item{bounds}{The bounds of the image, given as a four
  element array of \code{c(x_min, x_max, y_min, y_max)}. If
  bounds is NULL, it is calculated automatically from the
  centroids by giving a buffer region of \code{buffer}
  times the absolute spread of centroids.}

  \item{size}{Specify the x and y resolution of the image.
  Given as \code{c(nx, ny)}; defaults to c(500,500).}

  \item{buffer}{If \code{bounds} is NULL, then gives the
  height or width of the margins of the image containing
  the centroids.  For each x and y coordinates, this margin
  is equal to \code{buffer} times the difference between
  the minimum and maximum values present in the list of
  centroids.}
}
\value{
  A list with elements \code{stability}, \code{x},
  \code{y}, \code{bounds}, \code{centroids}, and
  \code{theta}. \code{stability} is the 2d image of size
  \code{size} to be plotted as the map of stable and
  unstable regions in the 2d space, \code{x} and \code{y}
  give the x and y positions in \code{stability},
  \code{theta} gives the original \code{theta} passed to
  the image, and \code{bounds} contains the \code{c(x_min,
  x_max, y_min, y_max)} bounds of the image.
}
\description{
  For a set of 2d centroids, shows the regions of stability
  and regions of instability for a given value of the
  perturbation hyperparameter.  The values in this plot
  indicate the contribution to the overall stability or
  instability from a point located at that value. This
  function is provided to demonstrate the intuitive
  behavior of the method and to help analyze 2d datasets.
}
\examples{
## Display the behavior of a set of centroids
library(easystab)

cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)

#to generate image with higher resolution, use larger size in the following line
Z <- make2dStabilityImage(cen, buffer=2, size=c(200,200))
image(Z$x, Z$y, Z$stability)
points(Z$centroids)

## Something more detailed; display how things change by theta

layout(matrix(1:4, ncol = 2, byrow=TRUE))
for(i in 1:4) {
  t <- (i - 1) * 0.5
  Z <- make2dStabilityImage(cen, theta=t, buffer=2, size=c(200,200))
  image(Z$x, Z$y, Z$stability, main = sprintf("Theta = \%1.2f.", t),
        xlab = "x", ylab="y")
}
}
\seealso{
  \code{\link{easystab}}
}

