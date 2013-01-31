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

\usage{
# Calculate perturbation stability of a list of clusterings
perturbationStability(clusterings, n_baselines = 32, seed = 0,
Kmap_mode = 0, theta = NULL)

# plot and summary of the stability report of a clustering as a heat map plot
plot(clustering, with_label = FALSE, classes = NULL, class_colors = NULL)
print(clustering)
summary(clustering)

# plot and summary of the stability scores produced by perturbationStability as a
  sequence of box plots
plot(clusterings)
print(clusterings)
summary(clusterings)

# plot pointwise stability image
plotStabilityImage(centroids, theta, image_nx, image_ny, image_x_lower =
0, image_y_upper = 0, image_y_lower = 0, image_y_upper =0)

}

\arguments{
  \item{clustering}{StabilityReport object. A clustering with stability analysis results
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
  \item{clusterings}{StabilitySequence object. The output of
    \code{perturbationStablity}
-- a list of clusters with perturbation stability analyses. }

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

## example with built in sample data
## clusterings in sample data already have 'dist' matrix included

library(easystab)

data(sample)
pl<-perturbationStability(clusterings)
res<-estimateK(pl)

l<-pl[[res$estimated_index]]

plot(pl)
plot(l)

## example with hclust function on iris data set

library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

dx <- dist(scale(X))
hc <- hclust(dx)
cl_list <- hclust_stability(dx, hc)

stability_sequence <- perturbationStability(cl_list)

print(stability_sequence)

summary(stability_sequence)

## example with kmeans function on iris data set

library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

km_list = list()

for(k in 1:10) {
  km_list[[k]] = kmeans(X, k)
}

stability_sequence <- perturbationStability(kmeans_stability(X, km_list))

# plots the sequence
plot(stability_sequence)

# plots the stability map
plot(stability_sequence[[3]], with_label = TRUE, classes = iris[,"Species"])

}

\seealso{
  \code{\link{perturbationStability}}, \code{\link{plotStabilityImage}}
}