\name{plotStabilityImage}

\alias{plotStabilityImage}

\title{Plots a picture of the relative regions of stability for a 2d clustering.}

\description{ For a set of 2d centroids, shows the regions of stability
  and regions of instability for a given value of the perturbation
  hyperparameter.  This function is provided to demonstrate the behavior
  of the method in an intuitive way.  }

\usage{
plotStabilityImage(centroids, theta = 0, bounds=NULL, size = c(501,501),
buffer=0.25)
}

\arguments{
  \item{centroids}{2D array of 2D centroid points, given as a K by 2
  array or matrix.}
  
  \item{theta}{The hyperparameter for the perturbation distribution of
  the prior.  This can be automatically computed by
  perturbationStability, or given manually.  Smaller values penalize the
  values on the boundary regions more severely.}

\item{bounds}{The bounds of the image, given as a four element array of
  \code{c(x_min, x_max, y_min, y_max)}. If bounds is NULL, it is
  calculated automatically from the centroids by giving a buffer region
  of \code{buffer} times the absolute spread of centroids.}

  \item{size}{Specify the x and y resolution of the image.  Given as
  \code{c(nx, ny)}.}
  
  \item{buffer}{If \code{bounds} is NULL, then gives the height or width
    of the margins of the image containing the centroids.  For each x
    and y coordinates, this margin is equal to \code{buffer} times the
    difference between the minimum and maximum values present in the
    list of centroids.}
}

  