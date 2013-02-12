StabilityColorMap <- colorRampPalette(c("black","red", "yellow", "white"))(512)


_processListOfClusterings <- function(clusterings) {
  ## Test whether a list of clusterings is actually a list of
  ## clusterings or a single clustering.

  if(! is.null(clusterings$dists)) {
    is_list <- FALSE

    clusterings$dists <- as.matrix(clusterings$dists)
    
    if(is.null(clusterings$labels)) {
      warning("'labels' attribute in clusterings not supplied; defaulting to minimum distance.")
      clusterings$labels <- apply(clusterings$dists, 1, which.min)
    }
  } else {
    is_list <- TRUE
    
    gave_warning <- FALSE
    
    for(i in 1:length(clusterings)) {
      if(is.null(clusterings[[i]]$dists)) {
        stop(sprintf("'dists' attribute not present in clustering component %d", i))
      }

      clusterings[[i]]$dists <- as.matrix(clusterings[[i]]$dists)
      
      if(is.null(clusterings[[i]]$labels)) {
        if(!gave_warning) {
          warning("'labels' attribute not supplied in all clusterings; defaulting to minimum distance.")
          gave_warning <- TRUE
        }
        
        clusterings$labels <- apply(clusterings[[i]]$dists, 1, which.min)
      }
    }
  }
  list(is_list <- is_list, clusterings <- clusterings)
}


f_theta <- function(t, clusterings, seed = 0, n_baselines = 32){
  score_total <- 0
  for(l in clusterings){
    X <- l$dists
    d <- dim(X)
    score_total <- score_total + .Call('_score', t(X), d[1], d[2], as.integer(seed), as.integer(n_baselines), t, FALSE, FALSE)
  }
  -score_total
}


#'Calculates the optimal prior parameter
#'
#'Calculate the optimal prior parameter theta by maximizing the difference in
#'overall stability aganst the baseline distributions. The theta parameter
#'indexes the strength of the perturbations, with smaller values translating
#'into stronger perturbations.
#'
#'
#'@param clusterings A single clustering or a list of clusterings. Each clustering of \code{n} data points into \code{K} clusters is specified primarily by matrix giving point to cluster distances.   Specifically, clustering must contain an of \code{n} by \code{K} distance matrix giving the point to cluster distance (\code{K} can be different across clusterings). Optionally, an array of \code{n} integer labels \code{labels} in \code{1,...,K} is expected; if not present, a warning is given and the labels are computed according to the minimum point-to-cluster distance.
#'
#'@param seed Random seed used for generating the baseline stability matrices.
#'
#'@param n_baselines The number of random baseline matrices to use in computing the stability scores.  Increase this number to get more accuracy at the expense of speed.
#'
#'@seealso \code{\link{easystab}}, \code{\link{perturbationStability}}


getOptTheta <- function(clusterings, seed = 0, n_baselines = 32){

  clusterings <- _processListOfClusterings(clusterings)$clusterings
  
  res <- optimize(f_theta, interval = c(-12, 12), tol = 0.00001,
                  clusterings = clusterings, seed = seed, n_baselines = n_baselines)
  
  res$minimum
}


#'Calculate clustering perturbation stability.
#'
#' Calculates the stability of clusterings under a non-parametric Bayesian perturbation as described in in [[PAPER]]. The exact method used is to perturb the cluster-to-point distances by scaling them with a shifted exponential random variable, then computing the probabilities of membership for each of the points under this perturbation.  This is compared against a set of randomly sampled bootstrapped baselines to determine the final stability score.
#' 
#'@param clusterings A single clustering or a list of clusterings. Each clustering of \code{n} data points into \code{K} clusters is specified primarily by matrix giving point to cluster distances.   Specifically, clustering must contain an of \code{n} by \code{K} distance matrix giving the point to cluster distance (\code{K} can be different across clusterings). Optionally, an array of \code{n} integer labels \code{labels} in \code{1,...,K} is expected; if not present, a warning is given and the labels are computed according to the minimum point-to-cluster distance.
#'
#'@param seed Random seed used for generating the baseline stability matrices.
#'
#'@param n_baselines The number of random baseline matrices to use in computing the stability scores.  Increase this number to get more accuracy at the expense of speed.
#'

#'@param Kmap_mode Whether to reorder the clusters in the stability map image for aesthetic reasons.  0 (default) means to not order them, 1 orders them by size, 2 orders them by average stability.  

#'@param theta The log of the rate parameter passed to the shifted exponential prior on the perturbations.  The parameter indexes the strength of the perturbations, with smaller values translating into stronger perturbations.  If NULL, theta is chosen by optimizing the overall stability against the baseline distributions as in \code{\link{getOptTheta}}.
                                        #'
#'@return Returns an object of type StabilitySequence if a list of clusterings is supplied, otherwise returns an object of type StabilityReport.  A StabilitySequence is essentially a list of StabilityReport objects corresponding to the original list of clusterings.
#'
#'A StabilityReport object contains the original \code{dists}, \code{labels} (possibly calculated), the scalar stability score \code{stability}, the empirical collection of stability scores \code{scores}, the theta parameter used or found \code{theta}, the individual membership probabilities of the points under perturbation, the \code{stability_matrix}
#'the sorted stability matrix used for plotting the behavior of the clustering.

perturbationStability <- function(clusterings, n_baselines = 32, seed = 0,
                                  Kmap_mode = 0, theta = NULL, test_pvalue = 0.05){
  if(seed < 0){
    warning("seed cannot be negative. your input is ", seed)
    return(NA)
  }
  
  is_list <- TRUE
  if(! is.null(names(clusterings))){
    tp_clustering <- clusterings
    clusterings <- list()
    clusterings[[1]] <- tp_clustering
    is_list <- FALSE
  }
  
  require(graphics)
  opt_theta <- theta
  
  if(is.null(theta)){
    opt_theta <- getOptTheta(clusterings, seed = seed, n_baselines = n_baselines)
  }

  n_points <- NULL

  for( idx in 1:length(clusterings)){
    l <- clusterings[[idx]]
    scores <- rep(as.numeric(NA), times = n_baselines)
    X <- l$dists
    
    if(is.null(X)){
      warning("clustering ", idx, " does not have point to centroid distance matrix\n")
      return(NA)
    }

    if(is.null(n_points))
      n_points <- nrow(X)
    else {
      if(n_points != nrow(X)) {
        warning("Inconsistent number of data points in list of clusterings!")
      }
    }

    d <- dim(X)
    .Call('_calculateScores', scores, t(X), d[1], d[2], as.integer(seed), as.integer(n_baselines), opt_theta, FALSE, FALSE)
    
    l$K <- d[2]
    l$stability <- mean(scores)
    l$stability_quantiles <- as.vector(quantile(scores, prob=c(0.025, 0.05, 0.95, 0.975), names=FALSE))
    l$scores <- scores
        
    Z <- -matrix(1, ncol = d[1], nrow = d[2])
    index_map <- rep(as.integer(NA), times=d[1])
    K_map <- rep(as.integer(NA), d[2])
    labels <- l$labels

    .Call('_sorted_stability_matrix', Z, index_map, K_map, t(X), labels, d[1], d[2], opt_theta, as.integer(Kmap_mode))

    l$sorted_stability_matrix <- t(Z)
    l$sorted_stability_matrix_index_map <- index_map
    l$sorted_stability_matrix_cluster_map <- K_map
    l$theta <- opt_theta

    class(l) <- "StabilityReport"

    clusterings[[idx]] <- l
  }

  class(clusterings) <- "StabilitySequence"
  
  if(is_list){
    ## Add in a bunch of estimates of the overall clustering stability

    cl_K_list = list()

    for(idx in 1:length(clusterings)){
      cl <- clusterings[[idx]]
      K <- cl$K
      
      cl$.original_index <- idx
      
      if(length(cl_K_list) < K
         ||is.null(cl_K_list[[K]])
         || cl_K_list[[K]]$stability < cl$stability) {
        
        cl_K_list[[K]] <- cl
        
      }
    }
      
      ## Get the most stable one.
    get_stability <- function(l) l$stability
    stability_vector <- sapply(cl_K_list, get_stability)
    e_idx <- which.max(stability_vector)
    
    ## See if there is not actually a most stable clustering. 
    if(cl_K_list[[e_idx]]$stability_quantiles[1] <= 0){

      ## There isn't a most stable clustering; 
      clusterings$estimated_K <- as.integer(1)
      clusterings$estimated_index <- as.integer(NA)

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
        
      clusterings$best <- l
      
    } else {
        
      estimated_K <- e_idx
        
      for(kidx in 1:(e_idx-1)) {
        if(!is.null(cl_K_list[[kidx]])) {
          res <- t.test(cl_K_list[[kidx]]$scores, cl_K_list[[e_idx]]$scores)
          t_stat <- res$statistic
          pv <- res$p.value * 0.5
          if(pv > test_pvalue){
            estimated_K <- kidx
            break
          }
        }
      }
      clusterings$estimated_K <- estimated_K
      clusterings$estimated_index <- cl_K_list[[estimated_K]]$.original_index
      clusterings$best <- cl_K_list[[estimated_K]]
    }
    return(clusterings)
  }else{
    return(clusterings[[1]])
  }
}



#'Adapts a single clustering, or list of clusterings, from \code{kmeans} to one
#'usable by \code{perturbationStability}.
#'
#'Given a list of clusterings, each from \code{kmeans}, returns a corresponding
#'list of clusterings suitable for input to \code{perturbationStability}.
#'
#'
#'@param x Matrix or data frame object containing the clustered data.
#'@param kmeans_output_list A list of clusterings, each being the output of the
#'kmeans function.
#'@return A list of clusterings that can be used as imput to the
#'\code{perturbationStability} function.
#'@seealso \code{\link{easystab}}, \code{\link{perturbationStability}},
#'\code{\link{hclust_stability}}
#'@examples
#'
#'## example with kmeans function on iris data set
#'
#'library(easystab)
#'
#'X <- iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]
#'
#'km_list = list()
#'
#'for(k in 1:10) {
#'  km_list[[k]] = kmeans(X, k)
#'}
#'
#'stability_sequence <- perturbationStability(kmeans_stability(X, km_list))
#'
#'# plots the sequence
#'plot(stability_sequence)
#'
#'# plots the stability map
#'plot(stability_sequence[[3]], with_label = TRUE, classes =
#'iris[,"Species"])
#'
#'############################################################
#'## Example with kmeans clustering on yeast data set
#'
#'yeast <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/yeast/yeast.data") 
#'
#'X <- yeast[,-c(1,10)]
#'class_label <- as.numeric(yeast[,10])
#'
#'# Sphere the data -- gives better results with 
#'X <- scale(X)
#'s <- svd(t(X))
#'X <- t(diag(1.0 / sqrt(s$d)) %*% t(s$u) %*% t(X))
#'
#'km_list = list() 
#'for(k in 1:12) {
#'  km_list[[k]] = kmeans(X, k, iter.max = 50, nstart=50)
#'}
#'
#'stability_sequence <- perturbationStability(kmeans_stability(X, km_list)) 
#'
#'print(stability_sequence)
#'
#'plot(stability_sequence)
#'
#'############################################################
#'## Example using kmeans_stability on a single clustering
#'
#'## Use X from previous yeast example
#'
#'## Works on a single clustering
#'km_cl <- kmeans(X, 8, iter.max = 50, nstart=50)
#'stability <- perturbationStability(kmeans_stability(X, km_cl))
#'
#'## Plot the stability -- a single clustering, so displays it as a
#'## stability map plot.
#'
#'plot(stability, with_label)
#'
#'
kmeans_stability <- function(x, kmeans_output_list) {

  is_list <- TRUE
  if(! is.null(names(kmeans_output_list))){
    tp_clustering <- kmeans_output_list
    kmeans_output_list <- list()
    kmeans_output_list[[1]] <- tp_clustering
    is_list <- FALSE
  }
  
  clusterings <- list()
  for(cl in kmeans_output_list) {
    cl$labels <- cl$cluster
    cl$dists <- rdist(x, cl$centers)
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
#'
#'@param dx Distance matrix as produced by \code{dists}, giving the
#'point-to-point distances.
#'@param hc Hierarchical clustering as produced by \code{hclust}.
#'@param clsnum_min Minimum cluster number (default 1).
#'@param clsnum_max Maximum cluster number (default 10).
#'@param method Method used to calculate the point-to-cluster distances from
#'the point-to-point distance matrix \code{dx} given.  Currently, the two
#'supported methods are "average", which takes the average of the distances
#'between the given point and all the points in the cluster (similar to average
#'linkage), and "median", which uses median distance to the points in the
#'cluster.
#'@return A list of clusterings suitable for use with
#'\code{perturbationStability}.
#'@seealso \code{\link{easystab}}, \code{\link{perturbationStability}},
#'\code{\link{clusterings_from_kmeans}}
hclust_stability <- function(dx, hc, clsnum_min = 1, clsnum_max = 10, method = "average"){

  if(method != "average" & method != "median"){
    warning("only average and median methods are supported")
    return(NA);
  }

  dx <- as.matrix(dx)
  
  n <- nrow(dx)
  
  if(n<clsnum_max){
    clsnum_max <- n
  }
  res <- cutree(hc, k=clsnum_min:clsnum_max)
  clusterings <- list()
  clsnum_min <- as.integer(clsnum_min)
  clsnum_max <- as.integer(clsnum_max)
  for(K in clsnum_min:clsnum_max){
    labels <- as.vector(res[,K-clsnum_min+1])
    dist_matrix <- matrix(as.numeric(NA), ncol = n, nrow = K)
    if(method == "average"){
      .Call('_calculateAverageLinkageDistances', dist_matrix, labels, n, K, dx)
    }else if (method == "median"){
      .Call('_calculateRepresentativeDistances', dist_matrix, labels, n, K, dx)
    }
    clusterings[[length(clusterings)+1]] <- list(dists = t(dist_matrix), labels = labels)
  }
  
  clusterings
}

print.StabilitySequence <- function(clusterings) {

  cat(sprintf("Perturbation Stability Sequence:\n"))
  cat(sprintf("  %d clusterings. \n", length(clusterings)))
  cat(sprintf("  Most stable clustering has %d clusters. \n", clusterings$estimated_K))
  
}

summary.StabilitySequence <- function(clusterings) {

  cat(sprintf("Perturbation Stability Sequence:\n"))
  cat(sprintf("  %d clusterings. \n", length(clusterings)))

  # Create a data frame to dump out the information
  
  X <- data.frame()
  
  idx <- 1
  for(idx in 1:length(clusterings)) {
    cl <- clusterings[[idx]]
    if(idx == clusterings$estimated_index)
      X[idx, " "] <- "BEST"
    else
      X[idx, " "] <- "    "
    
    X[idx, "K"] <- cl$K
    X[idx, "Score"] <- cl$stability
    X[idx, "95% CI"] <- sprintf("(%1.3f, %1.3f)",
                                cl$stability_quantiles[[1]],
                                cl$stability_quantiles[[4]])
  }
}



#'Plot the stability scores produced by perturbationStability as a sequence of
#'box plots.
#'
#'Summary display of the output of the perturbationStability function. Plot the
#'stability scores produced by perturbationStability as a sequence of box
#'plots.
#'
#'
#'@param clusterings The output of \code{perturbationStablity} -- a list of
#'clusters with perturbation stability analyses.
plot.StabilitySequence <- function(clusterings){
  get_scores <- function(l) l$scores
  get_K <- function(l) dim(l$dists)[2]
  score_list <- lapply(clusterings, get_scores)
  name_vector <- sapply(clusterings, get_K)
  names(score_list) <- name_vector
  boxplot(score_list)
}

print.StabilityReport <- function(clustering) {

  cat(sprintf("Perturbation Stability Report:\n"))
  cat(sprintf("  %d data points grouped into %d clusters.\n",
              dim(clustering$sorted_stability_matrix)[1], dim(clustering$K)))
  cat(sprintf("  Stability score = %f; 95%% confidence interval = [%f, %f].",
              clustering$stability,
              clustering$stability_quantiles[[1]],
              clustering$stability_quantiles[[4]]))
}

summary.StabilityReport <- function(clustering) {

  cat(sprintf("Perturbation Stability Summary:\n"))
  cat(sprintf("  %d data points grouped into %d clusters.\n",
              dim(clustering$sorted_stability_matrix)[1], dim(clustering$K)))
  cat(sprintf("  Stability score = %f; 95%% confidence interval = [%f, %f].",
              clustering$stability,
              clustering$stability_quantiles[[1]],
              clustering$stability_quantiles[[4]]))

  # Need more stuff in here; stability of each cluster for example...
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
#'
#'@param clustering A clustering with stability analysis results attached, as
#'given by an output of perturbationStability.
#'@param classes Auxiliary class labels for the data points, possibly from
#'known classes or other clusterings. The classes must integers in 1,...,L.  If
#'NULL, this column is not plotted.
#'@param class_colors Colors to use when plotting the auxiliary class labels.
#'If the given classes are in 1,...,L, it must be a list of at least L colors.
#'If NULL, \code{RColorBrewer} is used to choose representative colors.
#'Ignored if \code{classes} is \code{NULL}.
#'@seealso \code{\link{easystab}}
plot.StabilityReport <- function(clustering, classes = NULL, class_colors = NULL){
  
  require(grDevices)
  require(plotrix)
  
  sorted_stability_matrix <- clustering$sorted_stability_matrix

  nrow <- dim(sorted_stability_matrix)[1]
  ncol <- dim(sorted_stability_matrix)[2]

  color_func <- colorRampPalette(c("black","red", "yellow", "white"))
  color_map <-  color_func(512)

  mi <- min(sorted_stability_matrix)
  ma <- max(sorted_stability_matrix)
  
  matrix_colors <- apply(sorted_stability_matrix, 1:2, function(x) color_map[as.integer((x-mi)*511/(ma-mi))+1])

  if(is.null(classes)) {
    
    color2D.matplot(sorted_stability_matrix, cellcolors=matrix_colors, border=NA, xlab=NA, ylab=NA, axes=FALSE)
    axis(3, at=0.5:(ncol-0.5), labels=1:ncol)
    
  } else {
    labels <- as.integer(classes)

    index_map <- clustering$sorted_stability_matrix_index_map

    if(is.null(class_colors)){
      require(RColorBrewer)
      label_colors <- brewer.pal(max(labels), "Set3")
    }else{
      label_colors <- class_colors
    }
  
    label_column_colors <- rep(0, times = nrow)
  
    for(idx in 1:nrow){
      label_column_colors[idx] <- label_colors[labels[index_map[idx]]]
    }

    matrix_colors <- cbind(matrix_colors, label_column_colors)
    extended_matrix <- cbind(sorted_stability_matrix, labels)
    color2D.matplot(extended_matrix, cellcolors=matrix_colors, border=NA, xlab=NA, ylab=NA, axes=FALSE)
    axis(1, at=0.5:(ncol+0.5), labels = c(1:ncol, 'C'))

  }
}

makeStabilityImage <- function(centroids, theta = 0, bounds = NULL, size=c(501,501), buffer = 0.25) {

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
        t(centroids), nrow(centroids), theta, xvec, yvec, as.double(buffer))
  
  list(stability = stab_image,
       bounds = c(xvec[[1]], xvec[[size[[1]] ]], yvec[[1]], yvec[[size[[2]] ]]),
       size = size,
       centroids = centroids,
       theta = theta,
       x = xvec,
       y = yvec)
}



#'Plots a picture of the relative regions of stability for a 2d clustering.
#'
#'For a set of 2d centroids, shows the regions of stability and regions of
#'instability for a given value of the perturbation hyperparameter.  This
#'function is provided to demonstrate the behavior of the method in an
#'intuitive way.
#'
#'
#'@param centroids 2D array of 2D centroid points, given as a K by 2 array or
#'matrix.
#'@param theta The hyperparameter for the perturbation distribution of the
#'prior.  This can be automatically computed by perturbationStability, or given
#'manually.  Smaller values penalize the values on the boundary regions more
#'severely.
#'@param bounds The bounds of the image, given as a four element array of
#'\code{c(x_min, x_max, y_min, y_max)}. If bounds is NULL, it is calculated
#'automatically from the centroids by giving a buffer region of \code{buffer}
#'times the absolute spread of centroids.
#'@param size Specify the x and y resolution of the image.  Given as
#'\code{c(nx, ny)}.
#'@param buffer If \code{bounds} is NULL, then gives the height or width of the
#'margins of the image containing the centroids.  For each x and y coordinates,
#'this margin is equal to \code{buffer} times the difference between the
#'minimum and maximum values present in the list of centroids.
plotStabilityImage <- function(centroids, theta=0, bounds = NULL, size=c(501,501), buffer = 0.25) {

  res <- stab_image <- makeStabilityImage(centroids, theta, bounds, size, buffer)
  
  require(grDevices)
  require(plotrix)
  color2D.matplot(res$stab_image[nrow(res$stab_image):1,], border=NA, xlab=NA, ylab=NA, axes=FALSE)
  xlen <- length(res$xvec)
  ylen <- length(res$yvec)
  xarr <- c(1, as.integer((xlen+1)/4), as.integer((xlen+1)/2), as.integer((xlen+1)*3/4), xlen)
  yarr <- c(1, as.integer((ylen+1)/4), as.integer((ylen+1)/2), as.integer((ylen+1)*3/4), ylen)
  axis(1, at=xarr+0.5, labels=res$xvec[xarr])
  axis(2, at=yarr+0.5, labels=res$yvec[yarr])
}
