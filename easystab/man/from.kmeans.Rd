\name{from.kmeans}
\alias{from.kmeans}
\title{Adapts a single clustering, or list of clusterings, from \code{kmeans} to one usable by \code{perturbationStability}.}
\usage{
  from.kmeans(x, kmeans_output)
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
############################################################
## example with kmeans function on iris data set

library(easystab)

X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]

km_list = list()

for(k in 1:10) {
 km_list[[k]] = kmeans(X, k)
}

stability_collection <- perturbationStability(from.kmeans(X, km_list))

# plots the sequence
plot(stability_collection)

# plots the stability map of the 3 component case
plot(stability_collection[[3]], classes = iris[,"Species"])

############################################################
## Example with kmeans clustering on yeast data set

yeast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/yeast/yeast.data")

X <- yeast[,-c(1,10)]

# Sphere the data -- gives better results with the yeast data set
X <- scale(X)
s <- svd(t(X))
X <- t(diag(1.0 / sqrt(s$d)) \%*\% t(s$u) \%*\% t(X))

km_list = list()
for(k in 1:12) {
 km_list[[k]] = kmeans(X, k, iter.max = 50, nstart=50)
}

stability_collection <- perturbationStability(from.kmeans(X, km_list))

print(stability_collection)

plot(stability_collection)

############################################################
## Example using from.kmeans on a single clustering

## Use X from previous yeast example

## Works on a single clustering
km_cl <- kmeans(X, 8, iter.max = 50, nstart=50)
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

