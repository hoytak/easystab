\name{kmeans_stability}

\alias{kmeans_stability}

\title{Adapts a single clustering, or list of clusterings, from
  \code{kmeans} to one usable by \code{perturbationStability}.}

\description{Given a list of clusterings, each from \code{kmeans},
  returns a corresponding list of clusterings suitable for input to
  \code{perturbationStability}.}

\usage{ kmeans_stability(x, kmeans_output_list)}

\arguments{
  \item{x}{Matrix or data frame object containing the clustered data.}
  
  \item{kmeans_output_list}{A list of clusterings, each being the output of the
    kmeans function.}
}

\value{ A list of clusterings that can be used as imput to the
  \code{perturbationStability} function.}

\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}},
  \code{\link{hclust_stability}}
}

\examples{
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
plot(stability_sequence[[3]], with_label = TRUE, classes =
iris[,"Species"])

############################################################
## Example with kmeans clustering on yeast data set

yeast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/yeast/yeast.data") 

X <- yeast[,-c(1,10)]
class_label <- as.numeric(yeast[,10])

# Sphere the data -- gives better results with 
X <- scale(X)
s <- svd(t(X))
X <- t(diag(1.0 / sqrt(s$d)) %*% t(s$u) %*% t(X))

km_list = list() 
for(k in 1:12) {
  km_list[[k]] = kmeans(X, k, iter.max = 50, nstart=50)
}

stability_sequence <- perturbationStability(kmeans_stability(X, km_list)) 

print(stability_sequence)

plot(stability_sequence)

############################################################
## Example using kmeans_stability on a single clustering

## Use X from previous yeast example

## Works on a single clustering
km_cl <- kmeans(X, 8, iter.max = 50, nstart=50)
stability <- perturbationStability(kmeans_stability(X, km_cl))

## Plot the stability -- a single clustering, so displays it as a
## stability map plot.

plot(stability, with_label)

}