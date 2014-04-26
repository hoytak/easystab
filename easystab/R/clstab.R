#'@useDynLib easystab

StabilityColorMap <- function(n=512) { colorRampPalette(c("black","red", "yellow", "white"))(512)}

.processListOfClusterings <- function(clusterings) {
  ## Test whether a list of clusterings is actually a list of
  ## clusterings or a single clustering.

  give_label_warning <- TRUE

  # Function for checking the label stuff
  checked_labels <- function(clustering, suppress_warning, index = NULL) {

    X <- clustering$dists
    
    if(is.null(clustering$labels)) {
      if(!suppress_warning && give_label_warning) {
        warning("'labels' attribute in clustering not supplied; defaulting to minimum distance.")
        give_label_warning <<- FALSE
      }
      
      labels <- apply(X, 1, which.min)
      
    } else {
      labels <- as.integer(clustering$labels)

      if(max(labels) > ncol(X) || min(labels) < 1) {

        if(is.null(index))
          stop("'labels' vector contains indices out of range of 1,...,K.")
        else
          stop(sprintf("'labels' vector in entry %d contains indices out of range of 1,...,K.",
                       index))
      }
    }

    labels
  }

  if(is.matrix(clusterings)) {

    is_list <- FALSE

    X <- clusterings

    clusterings = list()
    clusterings$dists <- X

    n_points <- nrow(X)
    
    clusterings$labels <- checked_labels(clusterings, TRUE)

  } else if(! is.null(clusterings$dists)) {

    is_list <- FALSE

    clusterings$dists <- as.matrix(clusterings$dists)
    clusterings$labels <- checked_labels(clusterings, FALSE)

    n_points <- nrow(clusterings$dists)
    
  } else {

    is_list <- TRUE
    n_points <- NULL

    for(i in 1:length(clusterings)) {
      
      cl <- clusterings[[i]]
      
      if(is.matrix(cl)) {
        suppress_label_conversion_warning <- TRUE
        X <- cl
      } else {
        suppress_label_conversion_warning <- FALSE
      
        X <- clusterings[[i]]$dists
      
        if(is.null(X)) {
          stop("'dists' attribute not present in clustering component ", i, ".")
        }
      }

      if(is.null(n_points))
        n_points <- nrow(X)
      else {
        if(n_points != nrow(X)) {
          stop("Inconsistent number of data points in list of clusterings!")
        }
      }
      
      clusterings[[i]]$dists <- as.matrix(X)
      clusterings[[i]]$labels <- checked_labels(clusterings[[i]], suppress_label_conversion_warning)
    }
  }
  
  if(!is_list) {
    tp_clusterings <- clusterings
    clusterings <- list()
    clusterings[[1]] <- tp_clusterings
  }
  
  list(is_list = is_list,
       clusterings = clusterings,
       n_points = n_points)
}

f_theta <- function(theta, clusterings, seed, n_baselines){
  score_total <- 0
  for(idx in 1:length(clusterings)){
    l <- clusterings[[idx]]
    X <- l$dists
    labels <- l$labels
    d <- dim(X)
    score_total <- score_total + .Call('_score', t(X), labels, d[1], d[2],
                                       as.integer(seed), as.integer(n_baselines), theta, FALSE, FALSE)
  }
  
  print(c(theta, score_total))
  
  -score_total
}

## Orders the stability sequence by ensuring there is at most one for
## each K, and all in order
.orderedStabilityCollection <- function(clusterings, as_list_of_lists) {
  
  cl_K_list <- list()

  for(idx in 1:clusterings$n_clusterings){
    cl <- clusterings[[idx]]
    K <- cl$K
    
    if(as_list_of_lists) {
      if(length(cl_K_list) < K || is.null(cl_K_list[[K]])) {
        l <- list()
        l[[1]] <- cl
        cl_K_list[[K]] <- l
      } else {
        cl_K_list[[K]][[length(cl_K_list[[K]]) + 1]] <- cl
      }
    }
    else if(length(cl_K_list) < K
            ||  is.null(cl_K_list[[K]])
            || cl_K_list[[K]]$stability < cl$stability) {
      
      cl_K_list[[K]] <- cl
    }
  }
  
  cl_K_list <- Filter(Negate(is.null),  cl_K_list)

  if(as_list_of_lists) {
    sort_on_score <- function(cll) {
      scores <- as.vector(sapply(cll, function(cl) {cl$stability}))
      cll[order(scores,decreasing=TRUE)]
    }

    cl_K_list <- lapply(cl_K_list, sort_on_score)
  }

  cl_K_list
}

#'Calculates the optimal prior parameter
#'
#'Calculate the optimal prior parameter theta by maximizing the difference in overall stability aganst the baseline distributions. The theta parameter indexes the strength of the perturbations, with smaller values translating into stronger perturbations.
#'
#'
#'@param clusterings A single clustering or a list of clusterings. Each clustering of \code{n} data points into \code{K} clusters is specified primarily by matrix giving point to cluster distances.   Specifically, clustering must contain an of \code{n} by \code{K} distance matrix giving the point to cluster distance (\code{K} can be different across clusterings). Optionally, an array of \code{n} integer labels \code{labels} in \code{1,...,K} is expected; if not present, a warning is given and the labels are computed according to the minimum point-to-cluster distance.
#'
#'@param seed Random seed used for generating the baseline stability matrices.
#'
#'@param n_baselines The number of random baseline matrices to use in computing the stability scores.  Increase this number to get more accuracy at the expense of speed.
#'
#'@seealso \code{\link{easystab}}, \code{\link{perturbationStability}}
#'@export
#'
#'@examples
#'##################################################
#'## These examples produce exactly the same results as those in
#'## perturbationStability.
#'
#'library(easystab)
#'## Generate a fake dataset with 3 clusters
#'
#'cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
#'cl.size <- 100
#'X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]),
#'                   rnorm(cl.size,mean=cen[[1,2]])),
#'                         rbind(rnorm(cl.size,mean=cen[[2,1]]),
#'                               rnorm(cl.size,mean=cen[[2,2]])),
#'                         rbind(rnorm(cl.size,mean=cen[[3,1]]),
#'                               rnorm(cl.size,mean=cen[[3,2]]))))
#'dists  <- t(apply(X, 1, function(mu) {sqrt(rowSums((cen - mu)^2))}))
#'labels <- c(rep(1,100), rep(2,100), rep(3,100))
#'
#'## Takes same input as
#'theta <- getOptTheta(dists)
#'
#'## Apply to just the distance matrix
#'stability1 <- perturbationStability(dists, theta = theta)
#'
#'## Ways to display information
#'print(stability1)
#'summary(stability1)
#'plot(stability1, classes=labels)
#'
#'## Add in our labels
#'cl <- list(dists = dists, labels = labels)
#'stability2 <- perturbationStability(cl)
#'
#'print(stability2)
#'summary(stability2)
#'plot(stability2, classes=labels)
#'
#'## Now try several numbers of clusters using kmeans
#'km_list <- lapply(1:8, function(k) { kmeans(X, k, iter.max=20, nstart=30)})
#'cl_list <- from.kmeans(X, km_list)
#'stability_collection <- perturbationStability(cl_list)
#'
#'print(stability_collection)
#'summary(stability_collection)
#'plot(stability_collection)
#'
getOptTheta <- function(clusterings, seed = 0, n_baselines = 25){
  
  clusterings <- .processListOfClusterings(clusterings)$clusterings

  res <- optimize(f_theta, interval = c(0, 5), tol = 0.01 / n_baselines,
                  clusterings = clusterings, seed = seed,
                  n_baselines = n_baselines)
  
  res$minimum
}

#'Calculate clustering perturbation stability.
#'
#' Calculates the stability of clusterings under a non-parametric
#' Bayesian perturbation as described in in [[PAPER]]. The exact
#' method used is to perturb the cluster-to-point distances by
#' scaling them with a shifted exponential random variable, then
#' computing the probabilities of membership for each of the points
#' under this perturbation.  This is compared against a set of
#' randomly sampled bootstrapped baselines to determine the final
#' stability score.
#' 
#'@param clusterings A point-to-cluster distance matrix, a single
#'clustering, or a list of clusterings. Each clustering of \code{n}
#'data points into \code{K} clusters is specified primarily by
#'matrix giving point to cluster distances.  Specifically,
#'clustering must contain an of \code{n} by \code{K} distance
#'matrix giving the point to cluster distance (\code{K} can be
#'different across clusterings). Optionally, an array of \code{n}
#'integer labels \code{labels} in \code{1,...,K} is expected; if
#'not present, a warning is given and the labels are computed
#'according to the minimum point-to-cluster distance.  If only the
#'distance matrix is given, then points are assigned to their
#'minimum distance cluster with no warning.
#'
#'@param seed Random seed used for generating the baseline stability matrices.
#'
#'@param n_baselines The number of random baseline matrices to use
#'in computing the stability scores.  Increase this number to get
#'more accuracy at the expense of speed.
#'
#'@param theta The rate parameter passed to the shifted
#'exponential prior on the perturbations.  \code{theta} must
#'be non-negative; a warning is issued if \code{theta < 0}.
#'The parameter indexes
#'the strength of the perturbations, with smaller values
#'translating into stronger perturbations.  If NULL, theta is
#'chosen by optimizing the overall stability against the baseline
#'distributions as in \code{\link{getOptTheta}}.
#'
#'@param test_pvalue When selecting the best clustering among candidates
#'with a differing number of clusters, a one-sided t-test is performed
#'to choose the clustering having the smallest number of clusters and
#'statistically indistinguishable from the clustering with the highest
#'score. This is the level at which this t-test indicates that two stability
#'scores are statistically indistinguishable.
#'
#'@return Returns an object of type StabilityCollection if a list of
#'clusterings is supplied, otherwise returns an object of type
#'StabilityReport.  A StabilityCollection is essentially a list of
#'StabilityReport objects corresponding to the original list of
#'clusterings.
#'
#'A StabilityReport object contains the original \code{dists},
#'\code{labels} (possibly calculated), the scalar stability score
#'\code{stability}, the empirical collection of stability scores
#'\code{scores}, the theta parameter used or found \code{theta},
#'the individual membership probabilities of the points under
#'perturbation, the \code{stability_matrix} the sorted stability
#'matrix used for plotting the behavior of the clustering.
#'\code{print}, \code{summary}, and \code{plot} methods are
#'provided.
#'
#'@examples
#'## Generate a fake dataset with 3 clusters
#'cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
#'cl.size <- 100
## X is 300 x 2 matrix of 2d points, 100 from each of 3 components
#'X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
#'             rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
#'             rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))
#'dists  <- t(apply(X, 1, function(mu) {sqrt(rowSums((cen - mu)^2))}))
#'labels <- c(rep(1,100), rep(2,100), rep(3,100))
#'
#'## Apply to just the distance matrix
#'stability1 <- perturbationStability(dists)
#'
#'## Ways to display information
#'print(stability1)
#'summary(stability1)
#'plot(stability1, classes=labels)
#'
#'## Add in our labels
#'cl <- list(dists = dists, labels = labels)
#'stability2 <- perturbationStability(cl)
#'
#'print(stability2)
#'summary(stability2)
#'plot(stability2, classes=labels)
#'
#'## Now try several numbers of clusters using kmeans
#'km_list <- lapply(1:8, function(k) { kmeans(X, k, iter.max=20, nstart=30)})
#'cl_list <- from.kmeans(X, km_list)
#'stability_collection <- perturbationStability(cl_list)
#'
#'print(stability_collection)
#'summary(stability_collection)
#'plot(stability_collection)
#'
#'@seealso \code{\link{easystab}}, \code{\link{from.hclust}},
#'\code{\link{from.kmeans}}, \code{\link{getOptTheta}},
#'\code{\link{make2dStabilityImage}}
#'@export
perturbationStability <- function(clusterings, n_baselines = 25, seed = 0, theta = NULL, test_pvalue = 0.05){

  if(seed < 0){
    warning("Random seed cannot be negative. your input is ", seed)
    return(NA)
  }

  if(!is.null(theta) && theta < 0) {
    warning("Theta parameter cannot be negative; clipping to 0.")
    theta <- 0
  }

  cl_info <- .processListOfClusterings(clusterings)

  clusterings <- cl_info$clusterings

  is_list <- cl_info$is_list
  n_points <- cl_info$n_points

  if(is.null(theta)){
    res <- optimize(f_theta, interval = c(0, 5), tol = 0.01 / n_baselines,
                    clusterings = clusterings, seed = seed, n_baselines = n_baselines)
    opt_theta <- res$minimum
  } else {
    opt_theta <- theta
  }
  
  for( idx in 1:length(clusterings)){
    l <- clusterings[[idx]]
    scores <- rep(as.numeric(NA), times = n_baselines)
    X <- l$dists
    d <- dim(X)
    K <- d[[2]]
    
    confusion_matrix <- matrix(0, ncol = K, nrow = K)
    stability_matrix <- matrix(0, ncol = d[[1]], nrow = K)

    .Call('_calculateScores', scores, confusion_matrix, stability_matrix,
          t(X), l$labels, d[1], d[2],
          as.integer(seed), as.integer(n_baselines), opt_theta, FALSE, FALSE)
    
    l$K <- K
    l$stability <- mean(scores)
    l$stability_quantiles <- as.vector(quantile(scores, prob=c(0.025, 0.05, 0.95, 0.975), names=FALSE))
    l$scores <- scores

    l$confusion_matrix <- t(confusion_matrix)
    l$stability_matrix <- t(stability_matrix)
    l$theta <- opt_theta
    l$.original_index <- idx
    
    class(l) <- "StabilityReport"

    clusterings[[idx]] <- l
  }

  class(clusterings) <- "StabilityCollection"

  clusterings$n_clusterings = length(clusterings)

  if(is_list){
    ## Add in a bunch of estimates of the overall clustering stability
    cl_K_list <- .orderedStabilityCollection(clusterings, FALSE)
    
    ## Get the most stable one.
    stability_vector <- sapply(cl_K_list, function(l) { l$stability } )
    e_idx <- which.max(stability_vector)
    
    ## See if there is not actually a most stable clustering. 
    if(cl_K_list[[e_idx]]$stability_quantiles[1] <= 0){

      ## There isn't a most stable clustering; 
      clusterings$estimated_K <- as.integer(1)
      clusterings$best.index  <- NA

      ## Put a dummy clustering into the best slot
      l <- list()
      l$stability <- 0
      l$stability_quantiles <- as.vector(c(0,0,0,0))
      l$scores <- as.vector(rep(1, n_points))
      l$sorted_stability_matrix <-
          as.matrix(rep(1, n_points), nrow = n_points, ncol = 1)
      
      l$sorted_stability_matrix_index_map <- as.vector(1:n_points)
      l$sorted_stability_matrix_cluster_map <- as.vector(rep(1, n_points))
      l$theta <- 0

      class(l) <- "StabilityReport"
      
      clusterings$best <- l
      
    } else {
        
      estimated_K <- e_idx
        
      for(kidx in 1:(e_idx-1)) {
        res <- t.test(cl_K_list[[kidx]]$scores, cl_K_list[[e_idx]]$scores)
        t_stat <- res$statistic
        pv <- res$p.value * 0.5
        if(pv > test_pvalue){
          estimated_K <- kidx
          break
        }
      }
      
      clusterings$estimated_K <- estimated_K
      clusterings$best.index  <- cl_K_list[[estimated_K]]$.original_index
      clusterings$best        <- cl_K_list[[estimated_K]]
    }

    clusterings$theta <- opt_theta

    return(clusterings)
    
  } else {
    
    return(clusterings[[1]])
    
  }
}


#'Adapts a single clustering, or list of clusterings, from \code{kmeans} to one usable by \code{perturbationStability}.
#'
#'Given a clustering or list of clusterings, each from \code{kmeans}, returns a corresponding list of clusterings suitable for input to \code{perturbationStability}.
#'
#'@param X Matrix or data frame object containing the clustered data.  This is needed to compute the cluster to centroid distances.
#'
#'@param kmeans_output An output of kmeans objects, or list of such objects, each being the output of the
#'kmeans function.
#'
#'@return A clustering or list of clusterings that can be used as input to the
#'\code{perturbationStability} function.
#'
#'@seealso \code{\link{easystab}}, \code{\link{perturbationStability}},
#'\code{\link{from.hclust}}
#'
#'@examples
#'library(easystab)
#'
#'X <- scale(iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")])
#'
#'km_list <- lapply(1:12, function(k) { kmeans(X, k, iter.max=20, nstart=30)})
#'stability_collection <- perturbationStability(from.kmeans(X, km_list))
#'
#'## plots the sequence and stability map of the 3 component case
#'layout(matrix(1:2, nrow=1, ncol=2))
#'plot(stability_collection)
#'plot(stability_collection[[3]], classes = iris[,"Species"])
#'
#'############################################################
#'## Example with kmeans clustering on yeast data set
#'
#'yeast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/yeast/yeast.data")
#'
#'X <- scale(data.matrix(yeast[,-c(1,10)]))
#'
#'## To replicate results in paper, please comment out the following lines
#'## and increase the number of clusters considered to 12.  
#'rowmask <- yeast[,10] %in% c("MIT", "ME1", "ME2", "ME3")
#'yeast <- yeast[rowmask,]
#'X <- X[rowmask,]
#'
#'km_list <- lapply(1:6, function(k) { kmeans(X, k, iter.max=20, nstart=100)})
#'
#'stability_collection <- perturbationStability(from.kmeans(X, km_list))
#'
#'print(stability_collection)
#'
#'layout(matrix(1:2, nrow=1, ncol=2))
#'
#'## Plot the whole stability collection and stability map of the best one
#'plot(stability_collection)
#'plot(stability_collection$best, classes = yeast[,10])
#'
#'############################################################
#'## Example using from.kmeans on a single clustering
#'
#'## Use X from previous yeast example
#'
#'## Works on a single clustering
#'km_cl <- kmeans(X, 8, iter.max = 20, nstart=30)
#'stability <- perturbationStability(from.kmeans(X, km_cl))
#'
#'## Plot the stability -- a single clustering, so displays it as a
#'## stability map plot.
#'
#'plot(stability, classes=yeast[,10])
#'
#'@export
from.kmeans <- function(X, kmeans_output) {

  is_list <- TRUE
  if(! is.null(names(kmeans_output))){
    tp_clustering <- kmeans_output
    kmeans_output <- list()
    kmeans_output[[1]] <- tp_clustering
    is_list <- FALSE
  }
  
  clusterings <- list()
  for(cl in kmeans_output) {
    cl$labels <- cl$cluster
    cl$dists <- rdist(X, cl$centers)
    clusterings[[length(clusterings)+1]] <- cl
  }

  if(is_list){
    return(clusterings)
  }else{
    return(clusterings[[1]])
  }
}


#'Adapts the output of \code{hclust} for input into
#'\code{perturbationStability}.
#'
#'Adapts the output of \code{hclust} for use with \code{perturbationStability}
#'to give more information about the behavior of the hierarchical clustering
#'tree.
#'
#'@param dx Distance matrix as produced by \code{dists}, giving the
#'point-to-point distances.
#'
#'@param hc Hierarchical clustering as produced by \code{hclust}.
#'
#'@param k A list giving the numbers of clusters to cut the tree at;
#'this is passed to \code{\link{cutree}}.  Defaults to 1:10.
#'
#'@param method Method used to calculate the point-to-cluster distances from
#'the point-to-point distance matrix \code{dx} given.  Currently, the two
#'supported methods are "average", which takes the average of the distances
#'between the given point and all the points in the cluster (similar to average
#'linkage), and "median", which uses the median distance to the points in the
#'cluster.
#'
#'@return A list of clusterings suitable for use with
#'\code{perturbationStability}.
#'
#'@seealso \code{\link{easystab}}, \code{\link{perturbationStability}},
#'\code{\link{from.kmeans}}
#'
#'@examples
#'############################################################
#'## Interfacing with the hierarchical clustering method
#'library(easystab)
#'
#'## Generate a fake dataset with 3 clusters
#'cen <- matrix(c(0,-2,1,3,-3,1), ncol=2, byrow=TRUE)
#'cl.size <- 100
#'X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
#'            rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
#'            rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))
#'
#'dx <- dist(X)
#'hc <- hclust(dx)
#'cl_list <- from.hclust(dx,hc)
#'
#'stability_collection <- perturbationStability(cl_list)
#'
#'## Information about the stability sequence
#'print(stability_collection)
#'summary(stability_collection)
#'
#'## Plot the stability sequence
#'plot(stability_collection)
#'
#'############################################################
#'## A more detailed example using the UCI Wisconsin breast cancer dataset.
#'library(mlbench)
#'
#'# Load and cluster the Breast Cancer dataset using correlation distance.
#'data(BreastCancer)
#'
#'bcdata <- na.omit(BreastCancer)
#'
#'## Use 1 - (x %*% y) / (|x|_2 |y|_2) to compute divergence
#'X <- data.matrix(bcdata[,-c(1,11)])
#'Y <- X %*% t(X)
#'Ynorm <- diag(diag(Y)^(-1/2))
#'dx <- as.dist(1 - Ynorm %*% Y %*% Ynorm)
#'hc <- hclust(dx, method="complete")
#'
#'cl_list <- from.hclust(dx, hc, method = "median")
#'stability_collection <- perturbationStability(cl_list)
#'
#'# Information about the stability sequence
#'print(stability_collection)
#'summary(stability_collection)
#'
#'layout(matrix(1:2, nrow=1, ncol=2))
#'plot(stability_collection)
#'plot(stability_collection$best, classes = bcdata[,11])
#'@export
from.hclust <- function(dx, hc, k=1:10, method = "average") {

  if(method != "average" & method != "median"){
    warning("only average and median methods are supported")
    return(NA);
  }
  
  klist <- k
  dx <- as.vector(as.dist(dx))
  n <- as.integer((sqrt(1 + 8*length(dx)) + 1) / 2)
  res <- cutree(hc, k=klist)
  
  clusterings <- list()

  for(i in 1:length(klist)) {
    K <- klist[[i]]
    labels <- as.vector(res[,i])
    
    dist_matrix <- matrix(as.numeric(NA), ncol = n, nrow = K)
    if(method == "average"){
      .Call('_calculateAverageLinkageDistances', dist_matrix, labels, n, K, dx)
    }else if (method == "median"){
      .Call('_calculateRepresentativeDistances', dist_matrix, labels, n, K, dx)
    }
    clustering <- list()
    clustering$dists <- t(dist_matrix)
    clustering$labels <- labels
    clusterings[[length(clusterings)+1]] <- clustering
  }
  
  clusterings
}

#' Print a brief summary of the stability of a clustering collection.
#' 
#'@method print StabilityCollection
#'
#'@param x The output of \code{perturbationStablity} -- a
#'list of clusters with perturbation stability
#'analyses.
#'
#'@param ... optional arguments passed to internal functions
#'
#'@seealso \code{\link{easystab}}
#'@export
print.StabilityCollection <- function(x, ...) {

  clusterings = x
  
  cat(sprintf("Perturbation Stability Sequence:\n"))
  cat(sprintf("  %d clusterings. \n", clusterings$n_clusterings))
  cat(sprintf("  Most stable clustering has %d clusters. \n", clusterings$estimated_K))
}

#'Print a detaild summary of the stability of a clustering collection.
#'
#'@method summary StabilityCollection
#'
#'@param object The output of \code{perturbationStablity} -- a
#'list of clusters with perturbation stability
#'analyses.
#'
#'@param ... optional argements passed to internal functions
#'
#'@seealso \code{\link{easystab}}
#'@export
summary.StabilityCollection <- function(object, ...) {

  clusterings = object
  
  cat(sprintf("Perturbation Stability Sequence:\n"))
  cat(sprintf("  %d clusterings. \n", clusterings$n_clusterings))

  # Create a data frame to dump out the information
  
  X <- data.frame(matrix(0,nrow=clusterings$n_clusterings,ncol=4))
  colnames(X) <- c(" ", "K", "Score", "95% CI")

  
  
  for(idx in 1:clusterings$n_clusterings) {
    cl <- clusterings[[idx]]

    if(!is.na(clusterings$best.index) && idx == clusterings$best.index)
      X[[idx, " "]] <- "BEST"
    else
      X[[idx, " "]] <- "    "

    X[[idx, "K"]] <- cl$K
    X[[idx, "Score"]] <- cl$stability
    X[[idx, "95% CI"]] <- sprintf("(%1.3f, %1.3f)",
                                cl$stability_quantiles[[1]],
                                cl$stability_quantiles[[4]])
  }
  
  print(X)
}

#'Plot the stability scores produced by perturbationStability as a sequence of
#'box plots.
#'
#'Summary display of the output of the perturbationStability function. Plot the
#'stability scores produced by perturbationStability as a sequence of box
#'plots.
#'
#'@param x The output of \code{perturbationStablity} -- a
#'list of clusters with perturbation stability
#'analyses. Additionally:
#'
#'Set the \code{name} attribute of a specific clustering in order to
#'change the corresponding label on the box plots. For example,
#'\code{clusterings[[5]]$label <- "Clust5"} sets the displayed label
#'of that clustering, overriding the generated labels.
#'
#'Set the \code{color} attribute of a specific clustering in order to
#'change the color of boxplot.  The default is to color the "best"
#'one red, and the rest black (see \code{color.best} below).  This
#'overrides this behavior.
#'
#'@param sort Whether to sort the results in ascending order by the
#'number of clusters in the data, then by stability scores within the
#'clusters.  
#'
#'@param prune If sort is TRUE, and multiple clusterings are given
#'for a specific number of clusters, then show only the most stable
#'one from each group.  For example, if there were three clusterings
#'in the collection that had 5 clusters, only the most stable of
#'those three would be displayed.
#'
#'@param label.indices If \code{label.indices} is TRUE, then the
#'original indices from \code{clusterings} is included in the label
#'for each box plot; if FALSE, they are not included.  If
#'\code{label.indices} is NULL (default), then they are included only
#'if items in the graph are reordered.  Note that setting the
#'\code{label} attribute on the clusterings input overrides this.
#'
#'@param ... Additional parameters passed to the boxplot function.
#'See \code{\link{boxplot}} for more information.
#' 
#'@method plot StabilityCollection
#'@export
#'@examples
#'## Generate a fake dataset with 3 clusters
#'cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
#'cl.size <- 100
## X is 300 x 2 matrix of 2d points, 100 from each of 3 components
#'X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
#'             rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
#'             rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))
#'
#'
#'## Now try a range of numbers of clusters using kmeans
#'km_list1 <- lapply(1:6, function(k) { kmeans(X, k, iter.max=20, nstart=30)})
#'stabilities1 <- perturbationStability(from.kmeans(X, km_list1))
#'
#'plot(stabilities1)
#'
#'## Now plot each K with multiple runs of the clustering function.
#'## Now try several numbers of clusters using kmeans
#'km_list2 <- lapply(0:17, function(k) { kmeans(X, 1 + (k %% 6))})
#'stabilities2 <- perturbationStability(from.kmeans(X, km_list2))
#'
#'plot(stabilities2)
#'
#'## Plot the same thing, except without grouping by number of clusters
#'plot(stabilities2, sort=FALSE)
#'
#'## If two clusterings have the same number of clusters, plot only the
#'## most stable one.
#'plot(stabilities2, prune=TRUE, sort=FALSE)
#'
#'## Name the best one
#'stabilities2[[stabilities2$best.index]]$name <- "BEST!!!"
#'plot(stabilities2)
#'
#'@seealso \code{\link{easystab}}
plot.StabilityCollection <- function(x, sort = TRUE, prune = FALSE,
                                     label.indices = NULL, ...){

  clusterings = x
  n_cl <- clusterings$n_clusterings
  
  make.label <- function(cl, in_order) {
    if(!is.null(cl$name)) {
      return(cl$name)
    } else if(is.null(label.indices) && in_order) {
      s <- sprintf("%d", cl$K)
    } else if((is.null(label.indices) && !in_order) || label.indices) {
      s <- sprintf("%d/%d", cl$K, cl$.original_index)
    } else {
      s <- sprintf("%d", cl$K)
    }
    
    if(!is.na(clusterings$best.index)
       && cl$.original_index == clusterings$best.index) {
      return(sprintf("*%s*", s))
    } else {
      return(s)
    }
  }
  
  ylab <- "Stability Score"
  if(is.null(label.indices) || label.indices) {
    rotate_xlab <- TRUE
    xlab <- "# Clusters/Index"
  } else {
    rotate_xlab <- FALSE
    xlab <- "# Clusters"
  }

  if(sort && !prune) {

    cl_K_list <- .orderedStabilityCollection(clusterings, TRUE)

    ## TRUE if everything is in order, and the sequence of 
    if(length(cl_K_list) == clusterings$n_clusterings) {

      L1 <- lapply(cl_K_list, function(cl) { cl[[1]]$K } )
      L2 <- clusterings[[1]]$K:(clusterings[[clusterings$n_clusterings]]$K)
      
      scores <- lapply(cl_K_list, function(cll) { cll[[1]]$scores})

      if(length(L1) == length(L2) && all(L1 == L2)) {
        in_order <- TRUE
        xlab <- "# Clusters"
        rotate_xlab <- FALSE
      } else {
        in_order <- FALSE
      }
          
      names <- lapply(cl_K_list,
                      function(cll) { make.label(cll[[1]], in_order)})

      boxplot(scores, names = names, xlab = xlab, ylab = ylab,
              las=ifelse(rotate_xlab,3,0), outline=FALSE,...)
    
    } else {
      
      xpos <- 1
      xpos_vec <- NULL
      scores <- NULL
      names <- NULL
      
      for(cl_list in cl_K_list) {
        xpos_vec <- c(xpos_vec, xpos+(0:(length(cl_list) - 1)))
        xpos     <- xpos_vec[[length(xpos_vec)]] + 2
        scores   <- c(scores, lapply(cl_list, function(cl) {cl$scores}))
        names    <- c(names, lapply(cl_list, function(cl) {make.label(cl, FALSE)}))
      }
            
      boxplot(scores, names = names, xlab = xlab, ylab = ylab, 
              at=xpos_vec, las=ifelse(rotate_xlab,3,0),
              xlim=c(xpos_vec[[1]] - 1, xpos_vec[[length(xpos_vec)]]+1),
              outline=FALSE, 
              ...)
      
    }
  } else {
    
    if(sort && prune) 
      cl_list <- .orderedStabilityCollection(clusterings, FALSE)
    else
      cl_list <- lapply(1:n_cl, function(idx) {clusterings[[idx]]})
    
    scores <- lapply(cl_list, function(l) { l$scores} )
    name_vector <- lapply(cl_list, function(l) { make.label(l, TRUE)} )
    
    names(scores) <- name_vector

    boxplot(scores, main = "Adjusted Log-Stability Scores",
            xlab = "# Clusters", ylab = "Stability Score",
            las=ifelse(rotate_xlab,3,0), outline=FALSE,...)
  }

  yline(0)
}

#'Print a brief summary of the stability of an undividual clustering under perturbation. 
#'
#'@method print StabilityReport
#'
#'@param x A StabilityReport object, as given by an output
#'of perturbationStability.
#'
#'@param ... optional arguments passed to internal functions
#'
#'@seealso \code{\link{easystab}}
#'@export
print.StabilityReport <- function(x, ...) {

  clustering = x
  cat(sprintf("Perturbation Stability Report:\n"))
  cat(sprintf("  %d data points grouped into %d clusters.\n",
              nrow(clustering$dists), clustering$K))
  cat(sprintf("  Stability score = %f\n  95%% confidence interval = [%f, %f].\n",
              clustering$stability,
              clustering$stability_quantiles[[1]],
              clustering$stability_quantiles[[4]]))
}

#'Print a summary of the stability of an undividual clustering under perturbation. 
#'
#' Print a summary of the stability of an undividual clustering under
#' perturbation.  Summary includes individual cluster stabilites.
#'
#'@method summary StabilityReport
#'
#'@param object A StabilityReport object, as given by an output
#'of perturbationStability.
#'
#'@param ... optional arguments passed to internal functions
#'
#'@seealso \code{\link{easystab}}
#'@export
summary.StabilityReport <- function(object, ...) {

  clustering = object
  
  print.StabilityReport(clustering)

  cat(sprintf("Stability of Each Cluster Under Perturbation:\n\n"))

  K <- clustering$K
  C <- clustering$confusion_matrix
  X <- data.frame(C)

  colnames(X) <- lapply(1:K, function(k) {sprintf("--> %d", k)})

  Y <- data.frame(matrix(0,nrow=K,ncol=3))
  colnames(Y) <- c("Size", "Score", " ")
  
  for(idx in 1:K) {
    Y[[idx, "Size"]] <- sum(clustering$labels == idx)
    Y[[idx, "Score"]] <- C[[idx, idx]]
    Y[[idx, " "]] <- " "
  }

  print(cbind(Y, X), digits=3)
  cat(sprintf("\n  (Each row gives the average membership probabilties \n   under perturbation for data in that cluster.)\n\n"))
}

#'Display the stability of a clustering as a heat map plot.
#'
#'Plots the stability of a clustering as a heat map plot, showing the relative
#'stability of the different clusters, the data points, and the overall
#'behavior of the clustering.  The input is taken as a single clustering
#'analysis as given by perturbationStability.
#'
#'If \code{classes} are supplied (possibly from known classes or from another
#'clustering) version, they are plotted alongside the heatmap plot, with class
#'membership indexed by color.
#'
#'@param x A StabilityReport object, as given by an output
#'of perturbationStability.
#'
#'@param classes Auxiliary class labels for the data points, possibly from
#'known classes or other clusterings. The classes must be integers in 1,...,L.  If
#'NULL, this column is not plotted.
#'
#'@param class_colors Colors to use when plotting the auxiliary class labels.
#'If the given classes are in 1,...,L, it must be a list of at least L colors.
#'If NULL, \code{RColorBrewer} is used to choose representative colors.
#'Ignored if \code{classes} is \code{NULL}.
#'
#'@param sort.clusters Whether to sort the clusters in the stability
#'map image for aesthetic reasons.  0 (default) means to not reorder
#'them, 1 orders them by cluster size, and 2 orders them by average stability.
#'
#'@param ... optional arguments passed to internal functions
#'
#'@seealso \code{\link{easystab}}
#'
#'@method plot StabilityReport
#'@export
#'@examples
#'## Generate a fake dataset with 3 clusters
#'cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
#'cl.size <- 100
## X is 300 x 2 matrix of 2d points, 100 from each of 3 components
#'X <- t(cbind(rbind(rnorm(cl.size,mean=cen[[1,1]]), rnorm(cl.size,mean=cen[[1,2]])),
#'             rbind(rnorm(cl.size,mean=cen[[2,1]]), rnorm(cl.size,mean=cen[[2,2]])),
#'             rbind(rnorm(cl.size,mean=cen[[3,1]]), rnorm(cl.size,mean=cen[[3,2]]))))
#'dists  <- t(apply(X, 1, function(mu) {sqrt(rowSums((cen - mu)^2))}))
#'labels <- c(rep(1,100), rep(2,100), rep(3,100))
#'
#'## Apply to just the distance matrix
#'stability1 <- perturbationStability(dists)
#'
#'## Ways to display information
#'print(stability1)
#'summary(stability1)
#'plot(stability1, classes=labels)
#'
#'## Add in our labels
#'cl <- list(dists = dists, labels = labels)
#'stability2 <- perturbationStability(cl)
#'
#'print(stability2)
#'summary(stability2)
#'plot(stability2, classes=labels)
#'
plot.StabilityReport <- function(x, classes = NULL, class_colors = NULL, sort.clusters = 0, ...){

  clustering = x
  
  require(grDevices)
  require(plotrix)

  n <- nrow(clustering$stability_matrix)
  K <- ncol(clustering$stability_matrix)
  sorted_stability_matrix <- matrix(-1, ncol = n, nrow = K)
  index_map <- rep(as.integer(NA), n)
  K_map <- rep(as.integer(NA), K)
  labels <- clustering$labels

  .Call('_sort_stability_matrix', sorted_stability_matrix, index_map, K_map,
        t(clustering$stability_matrix), labels,
        n, K, as.integer(sort.clusters))

  sorted_stability_matrix <- t(sorted_stability_matrix)
  
  color_func <- colorRampPalette(c("black","red", "yellow", "white"))
  color_map <-  color_func(512)
  
  matrix_colors <- apply(sorted_stability_matrix, 1:2, function(x) color_map[as.integer(x*511)+1])

  if(is.null(classes)) {
    
    color2D.matplot(sorted_stability_matrix, cellcolors=matrix_colors, border=NA, xlab=NA, ylab=NA, axes=FALSE)
    axis(3, at=0.5:(n-0.5), labels=1:n)
    
  } else {
    classes <- as.integer(classes)
    n_classes <- max(classes)

    if(length(classes) != n)
      stop(sprintf("Length of supplied classes (%d) does not match the number of data points (%d).",
                   length(classes), n))
    
    if(is.null(class_colors)){
      require(RColorBrewer)
      label_colors <- brewer.pal(max(c(n_classes, 3)), "Set3")
    }else{
      label_colors <- class_colors
    }
  
    label_column_colors <- rep(0, n_classes)

    for(idx in 1:n) {
      label_column_colors[[idx]] <- label_colors[[classes[index_map[idx]]]]
    }

    matrix_colors <- cbind(matrix_colors, label_column_colors)
    extended_matrix <- cbind(sorted_stability_matrix, classes)
    color2D.matplot(extended_matrix, cellcolors=matrix_colors, border=NA, xlab=NA, ylab=NA, axes=FALSE)
    axis(1, at=0.5:(K+0.5), labels = c(1:K, 'C'))
  }
}

#'Creates an image of the relative regions of stability for a 2d clustering.
#'
#'For a set of 2d centroids, shows the regions of stability and
#'regions of instability for a given value of the perturbation
#'hyperparameter.  The values in this plot indicate the
#'contribution to the overall stability or instability from a point
#'located at that value. This function is provided to demonstrate
#'the intuitive behavior of the method and to help analyze 2d
#'datasets.
#'
#'@param centroids Array of 2D centroid points, given as a K by 2 array or
#'matrix.
#'
#'@param theta The rate parameter passed to the shifted
#'exponential prior on the perturbations.  \code{theta} must
#'be non-negative; a warning is issued if \code{theta < 0}.
#'The parameter indexes
#'the strength of the perturbations, with smaller values
#'translating into stronger perturbations.  If NULL, theta is
#'chosen by optimizing the overall stability against the baseline
#'distributions as in \code{\link{getOptTheta}}.
#'
#'@param bounds The bounds of the image, given as a four element array of
#'\code{c(x_min, x_max, y_min, y_max)}. If bounds is NULL, it is calculated
#'automatically from the centroids by giving a buffer region of \code{buffer}
#'times the absolute spread of centroids.
#'
#'@param size Specify the x and y resolution of the image.  Given as
#'\code{c(nx, ny)}; defaults to c(500,500).
#'
#'@param buffer If \code{bounds} is NULL, then gives the height or width of the
#'margins of the image containing the centroids.  For each x and y coordinates,
#'this margin is equal to \code{buffer} times the difference between the
#'minimum and maximum values present in the list of centroids.
#'
#'@return A list with elements \code{stability}, \code{x},
#'\code{y}, \code{bounds}, \code{centroids}, and \code{theta}.
#'\code{stability} is the 2d image of size \code{size} to be
#'plotted as the map of stable and unstable regions in the 2d
#'space, \code{x} and \code{y} give the x and y positions in
#'\code{stability}, \code{theta} gives the original \code{theta}
#'passed to the image, and \code{bounds} contains the
#'\code{c(x_min, x_max, y_min, y_max)} bounds of the image.
#'
#'@examples
#'## Display the behavior of a set of centroids 
#'library(easystab)
#'
#'cen <- matrix(c(0,-2,1,2,-2,1), ncol=2, byrow=TRUE)
#'
#'#to generate image with higher resolution, use larger size in the following line
#'Z <- make2dStabilityImage(cen, buffer=2, size=c(200,200))
#'image(Z$x, Z$y, Z$stability)
#'points(Z$centroids)
#'
#'## Something more detailed; display how things change by theta
#' 
#' layout(matrix(1:4, ncol = 2, byrow=TRUE))
#' for(i in 1:4) {
#'   t <- (i - 1) * 0.5
#'   Z <- make2dStabilityImage(cen, theta=t, buffer=2, size=c(200,200))
#'   image(Z$x, Z$y, Z$stability, main = sprintf("Theta = %1.2f.", t),
#'         xlab = "x", ylab="y")
#' }
#' 
#'@export
#'
#'@seealso \code{\link{easystab}}
make2dStabilityImage <- function(centroids, theta = 1, bounds = NULL, size=c(500,500), buffer = 0.25) {

  if(!is.null(theta) && theta < 0) {
    warning("Theta parameter cannot be negative; clipping to 0.")
    theta <- 0
  }
  
  image_nx <- size[[1]]
  image_ny <- size[[2]]

  if(is.null(bounds)) {
    image_x_lower <- 0
    image_x_upper <- 0
    image_y_lower <- 0
    image_y_upper <- 0
  } else {
    image_x_lower <- bounds[[1]]
    image_x_upper <- bounds[[2]]
    image_y_lower <- bounds[[3]]
    image_y_upper <- bounds[[4]]
  } 
  
  stab_image <- matrix(as.numeric(NA), ncol = image_ny, nrow = image_nx)
  xvec <- rep(as.numeric(NA), times = as.integer(image_nx))
  yvec <- rep(as.numeric(NA), times = as.integer(image_ny))
  .Call('_make_stability_image', stab_image, as.integer(image_nx), as.integer(image_ny),
        image_x_lower, image_x_upper, image_y_lower, image_y_upper,
        t(centroids), nrow(centroids), as.double(theta), xvec, yvec, as.double(buffer))
  
  list(stability = stab_image,
       bounds = c(xvec[[1]], xvec[[size[[1]] ]], yvec[[1]], yvec[[size[[2]] ]]),
       size = size,
       centroids = centroids,
       theta = theta,
       x = xvec,
       y = yvec)
}
