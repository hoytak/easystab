\name{plot.StabilityCollection}
\alias{plot.StabilityCollection}
\title{Plot the stability scores produced by perturbationStability as a sequence of
box plots.}
\usage{
  \method{plot}{StabilityCollection} (x, sort = TRUE,
    prune = FALSE, label.indices = NULL, ...)
}
\arguments{
  \item{x}{The output of \code{perturbationStablity} -- a
  list of clusters with perturbation stability analyses.
  Additionally:

  Set the \code{name} attribute of a specific clustering in
  order to change the corresponding label on the box plots.
  For example, \code{clusterings[[5]]$label <- "Clust5"}
  sets the displayed label of that clustering, overriding
  the generated labels.

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
}
\description{
  Summary display of the output of the
  perturbationStability function. Plot the stability scores
  produced by perturbationStability as a sequence of box
  plots.
}
\examples{
## Generate a fake dataset with 3 clusters
cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
cl.size <- 100
X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
            rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
            rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))


## Now try a range of numbers of clusters using kmeans
km_list1 <- lapply(1:6, function(k) { kmeans(X, k, iter.max=20, nstart=30)})
stabilities1 <- perturbationStability(from.kmeans(X, km_list1))

plot(stabilities1)

## Now plot each K with multiple runs of the clustering function.
## Now try several numbers of clusters using kmeans
km_list2 <- lapply(0:17, function(k) { kmeans(X, 1 + (k \%\% 6))})
stabilities2 <- perturbationStability(from.kmeans(X, km_list2))

plot(stabilities2)

## Plot the same thing, except without grouping by number of clusters
plot(stabilities2, sort=FALSE)

## If two clusterings have the same number of clusters, plot only the
## most stable one.
plot(stabilities2, prune=TRUE, sort=FALSE)

## Name the best one
stabilities2[[stabilities2$best.index]]$name <- "BEST!!!"
plot(stabilities2)
}
\seealso{
  \code{\link{easystab}}
}

