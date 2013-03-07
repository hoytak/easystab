\name{from.kmeans}
\alias{from.kmeans}
\title{Adapts a single clustering, or list of clusterings, from \code{kmeans} to one usable by \code{perturbationStability}.}
\usage{
  from.kmeans(X, kmeans_output)
}
\arguments{
  \item{X}{Matrix or data frame object containing the
  clustered data.  This is needed to compute the cluster to
  centroid distances.}

  \item{kmeans_output}{An output of kmeans objects, or list
  of such objects, each being the output of the kmeans
  function.}
}
\value{
  A clustering or list of clusterings that can be used as
  input to the \code{perturbationStability} function.
}
\description{
  Given a clustering or list of clusterings, each from
  \code{kmeans}, returns a corresponding list of
  clusterings suitable for input to
  \code{perturbationStability}.
}
\examples{
library(easystab)

X <- scale(iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")])

km_list <- lapply(1:12, function(k) { kmeans(X, k, iter.max=20, nstart=30)})
stability_collection <- perturbationStability(from.kmeans(X, km_list))

## plots the sequence and stability map of the 3 component case
layout(matrix(1:2, nrow=1, ncol=2))
plot(stability_collection)
plot(stability_collection[[3]], classes = iris[,"Species"])

############################################################
## Example with kmeans clustering on yeast data set

yeast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/yeast/yeast.data")

X <- scale(data.matrix(yeast[,-c(1,10)]))

## To replicate results in paper, please comment out the following lines
## and increase the number of clusters considered to 12.
rowmask <- yeast[,10] \%in\% c("MIT", "ME1", "ME2", "ME3")
yeast <- yeast[rowmask,]
X <- X[rowmask,]

km_list <- lapply(1:6, function(k) { kmeans(X, k, iter.max=20, nstart=100)})

stability_collection <- perturbationStability(from.kmeans(X, km_list))

print(stability_collection)

layout(matrix(1:2, nrow=1, ncol=2))

## Plot the whole stability collection and stability map of the best one
plot(stability_collection)
plot(stability_collection$best, classes = yeast[,10])

############################################################
## Example using from.kmeans on a single clustering

## Use X from previous yeast example

## Works on a single clustering
km_cl <- kmeans(X, 8, iter.max = 20, nstart=30)
stability <- perturbationStability(from.kmeans(X, km_cl))

## Plot the stability -- a single clustering, so displays it as a
## stability map plot.

plot(stability, classes=yeast[,10])
}
\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}},
  \code{\link{from.hclust}}
}

