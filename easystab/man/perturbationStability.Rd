\name{perturbationStability}
\alias{perturbationStability}
\title{Calculate clustering perturbation stability.}
\usage{
  perturbationStability(clusterings, n_baselines = 32,
    seed = 0, theta = NULL, test_pvalue = 0.05)
}
\arguments{
  \item{clusterings}{A point-to-cluster distance matrix, a
  single clustering, or a list of clusterings. Each
  clustering of \code{n} data points into \code{K} clusters
  is specified primarily by matrix giving point to cluster
  distances.  Specifically, clustering must contain an of
  \code{n} by \code{K} distance matrix giving the point to
  cluster distance (\code{K} can be different across
  clusterings). Optionally, an array of \code{n} integer
  labels \code{labels} in \code{1,...,K} is expected; if
  not present, a warning is given and the labels are
  computed according to the minimum point-to-cluster
  distance.  If only the distance matrix is given, then
  points are assigned to their minimum distance cluster
  with no warning.}

  \item{seed}{Random seed used for generating the baseline
  stability matrices.}

  \item{n_baselines}{The number of random baseline matrices
  to use in computing the stability scores.  Increase this
  number to get more accuracy at the expense of speed.}

  \item{theta}{The log of the rate parameter passed to the
  shifted exponential prior on the perturbations.  The
  parameter indexes the strength of the perturbations, with
  smaller values translating into stronger perturbations.
  If NULL, theta is chosen by optimizing the overall
  stability against the baseline distributions as in
  \code{\link{getOptTheta}}.}
}
\value{
  Returns an object of type StabilityCollection if a list
  of clusterings is supplied, otherwise returns an object
  of type StabilityReport.  A StabilityCollection is
  essentially a list of StabilityReport objects
  corresponding to the original list of clusterings.

  A StabilityReport object contains the original
  \code{dists}, \code{labels} (possibly calculated), the
  scalar stability score \code{stability}, the empirical
  collection of stability scores \code{scores}, the theta
  parameter used or found \code{theta}, the individual
  membership probabilities of the points under
  perturbation, the \code{stability_matrix} the sorted
  stability matrix used for plotting the behavior of the
  clustering. \code{print}, \code{summary}, and \code{plot}
  methods are provided.
}
\description{
  Calculates the stability of clusterings under a
  non-parametric Bayesian perturbation as described in in
  [[PAPER]]. The exact method used is to perturb the
  cluster-to-point distances by scaling them with a shifted
  exponential random variable, then computing the
  probabilities of membership for each of the points under
  this perturbation.  This is compared against a set of
  randomly sampled bootstrapped baselines to determine the
  final stability score.
}
\examples{
## Toy example using a simple mixture model
cen <- matrix(c(0,0,1,2,-1,0.5), ncol=2)

## X is 300 x 2 matrix of 2d points, 100 from each of 3 components
X <- t(cbind(rbind(rnorm(100,mean=cen[[1,1]],sd=1),rnorm(100,mean=cen[[1,2]],sd=1)),
            rbind(rnorm(100,mean=cen[[2,1]],sd=2),rnorm(100,mean=cen[[2,2]],sd=2)),
            rbind(rnorm(100,mean=cen[[3,1]],sd=1),rnorm(100,mean=cen[[3,2]],sd=1))))

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

## Now try several numbers of clusters using kmeans
km_list = list()

for(k in 1:8)
 km_list[[k]] = kmeans(X, k, iter.max=50, nstart=100)

cl_list <- from.kmeans(X, km_list)
stability_sequence <- perturbationStability(cl_list)

print(stability_sequence)
summary(stability_sequence)
plot(stability_sequence)
}

