f_theta <- function(t, clusterings, seed = 0, n_baselines = 32){
  score_total <- 0
  for(l in clusterings){
    X <- l$dists
    d <- dim(X)
    score_total <- score_total + .Call('_score', t(X), d[1], d[2], as.integer(seed), as.integer(n_baselines), t, FALSE, FALSE)
  }
  -score_total
}

getOptTheta <- function(clusterings, seed = 0, n_baselines = 32){
    res <- optimize(f_theta, interval = c(-12, 12), tol = 0.00001, clusterings = clusterings, seed = seed, n_baselines = n_baselines)
    res$minimum
}

make_stability_image <- function(centroids, theta, image_nx, image_ny, image_x_lower = 0, image_x_upper = 0, image_y_lower = 0, image_y_upper = 0){
  stab_image <- matrix(as.numeric(NA), ncol = image_ny, nrow = image_nx)
  xvec <- rep(as.numeric(NA), times = as.integer(image_nx))
  yvec <- rep(as.numeric(NA), times = as.integer(image_ny))
  .Call('_make_stability_image', stab_image, as.integer(image_nx), as.integer(image_ny), image_x_lower, image_x_upper, image_y_lower, image_y_upper, t(centroids), nrow(centroids), theta, xvec, yvec)
  list(stab_image = t(stab_image), xvec = xvec, yvec = yvec)
}

perturbationStability <- function(clusterings, n_baselines = 32, seed = 0, Kmap_mode = 0, theta = NULL){
  require(graphics)
  opt_theta = theta
  if(is.null(theta)){
#    res <- optimize(f_theta, interval = c(-8, 8), tol = 0.00001, clusterings = clusterings, n_baselines = n_baselines)
#    opt_theta <- res$minimum
    opt_theta <- getOptTheta(clusterings, seed = seed, n_baselines = n_baselines)
  }
  for( idx in 1:length(clusterings)){
    l <- clusterings[[idx]]
    scores = rep(as.numeric(NA), times = n_baselines)
    X <- l$dists
    d <- dim(X)
    .Call('_calculateScores', scores, t(X), d[1], d[2], as.integer(seed), as.integer(n_baselines), opt_theta, FALSE, FALSE)
    l$stability <- mean(scores)
    l$stability_quantiles <- as.vector(quantile(scores, prob=c(0.025, 0.05, 0.95, 0.975), names=FALSE))
    l$scores = scores
    
    Z <- -matrix(1, ncol = d[1], nrow = d[2])
    index_map <- rep(as.integer(NA), times=d[1])
    K_map <- rep(as.integer(NA), d[2])
    labels <- l$labels

    .Call('_sorted_stability_matrix', Z, index_map, K_map, t(X), labels, d[1], d[2], opt_theta, as.integer(Kmap_mode))

    l$sorted_stability_matrix <- t(Z)
    l$sorted_stability_matrix_index_map <- index_map
    l$sorted_stability_matrix_cluster_map <- K_map
    l$theta <- opt_theta
    clusterings[[idx]] <- l
  }

  clusterings
}

gen_kmeans_clusterings <- function(clustering_size = 9){
  require(fields)
  x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
        matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
  clusterings = list()
  for(i in 2:(clustering_size+1)){
     cl <- kmeans(x, i)
     cl$labels = cl$cluster
     cl$dists <- rdist(x, cl$centers)
     clusterings[[length(clusterings)+1]] = cl
  }
  clusterings
}

kmeans_stability <- function(x, clusterings) {

  for(cl in clusterings) {
    cl$labels = cl$cluster
    cl$dists <- rdist(x, cl$centers)
  }
  
  clusterings
}

hclust_stability <- function(dx, hc, clsnum_min = 2, clsnum_max = 10, method = "average"){
  n <- nrow(dx)
  
  if(n<clsnum_max){
    clsnum_max <- n
  }
  res <- cutree(hc, k=clsnum_min:clsnum_max)
  clusterings = list()
  clsnum_min = as.integer(clsnum_min)
  clsnum_max = as.integer(clsnum_max)
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

estimateK <- function(clusterings, p_value = 0.05){
  get_stability <- function(l) l$stability
  stability_vector <- sapply(clusterings, get_stability)
  
  e_idx = which.max(stability_vector)
  for(idx in 1:(e_idx-1)){
    res <- t.test(clusterings[[idx]]$scores, clusterings[[e_idx]]$scores)
    t_stat <- res$statistic
    pv <- res$p.value * 0.5
    if(pv > p_value){
      break
    }
  }
  if(clusterings[[e_idx]]$stability_quantiles[1] <=0 ){
    estimated_K = as.integer(1)
    estimated_index = as.integer(NA)
  }else{
   estimated_K = dim(clusterings[[e_idx]]$dists)[2]
   estimated_index = e_idx
  }
  list(estimated_K = estimated_K, estimated_index = estimated_index)
}

plotStabilitySequence <- function(clusterings){
  get_scores <- function(l) l$scores
  get_K <- function(l) dim(l$dists)[2]
  score_list <- lapply(clusterings, get_scores)
  name_vector <- sapply(clusterings, get_K)
  names(score_list) <- name_vector
  boxplot(score_list)
}

#plotStabilityMapO <- function(clustering){
#  sorted_stability_matrix = clustering$sorted_stability_matrix
#  require(plotrix)
#  cellcol<-color.scale(sorted_stability_matrix,c(0,1),0,0)
#  color2D.matplot(sorted_stability_matrix, cellcolors=cellcol, xlab=NA, ylab=NA, axes=FALSE)
#}

plotStabilityMap <- function(clustering){
  sorted_stability_matrix = clustering$sorted_stability_matrix
  ncol = dim(sorted_stability_matrix)[2]
  require(grDevices)
  require(plotrix)
  colors <- colorRampPalette(c("black","red", "yellow", "white"))
  mi = 0 #min(sorted_stability_matrix)
  ma = 1 #max(sorted_stability_matrix)
  mx <- apply(sorted_stability_matrix, 1:2, colors)
  color2D.matplot(sorted_stability_matrix, cellcolors=mx, border=NA, xlab=NA, ylab=NA, axes=FALSE)
  axis(3, at=0.5:(ncol-0.5), labels=1:ncol)
}

plotLabeledStabilityMap <- function(clustering, labels, label_colors = NULL){
  require(grDevices)
  require(plotrix)
  
  sorted_stability_matrix <- clustering$sorted_stability_matrix
  index_map <- clustering$sorted_stability_matrix_index_map
  
  label_colors <- brewer.pal(max(labels), "Set3")

  if(is.null(label_colors)){
    require(RColorBrewer)
    colors <- colorRampPalette(c("black","red", "yellow", "white"))
  }else{
    colors <- label_colors
  }
  
  mi = min(sorted_stability_matrix)
  ma = max(sorted_stability_matrix)
  matrix_colors <- apply(sorted_stability_matrix, 1:2, colors)

  nrow = dim(sorted_stability_matrix)[1]
  ncol = dim(sorted_stability_matrix)[2]
  
  label_column_colors= rep(0, times = nrow)
  
  for(idx in 1:nrow){
    label_column_colors[idx] <- label_colors[labels[index_map[idx]]]
  }

  matrix_colors <- cbind(matrix_colors, label_column_colors)
  extended_matrix <- cbind(sorted_stability_matrix, labels)
  color2D.matplot(extended_matrix, cellcolors=matrix_colors, border=NA, xlab=NA, ylab=NA, axes=FALSE)
  axis(1, at=0.5:(ncol+0.5), labels=c(1:ncol, 'C'))

  colors
}

plotStabilityImage <- function(centroids, theta, image_nx, image_ny, image_x_lower = 0, image_x_upper = 0, image_y_lower = 0, image_y_upper = 0){
  res <- stab_image <- make_stability_image(centroids, theta, image_nx, image_ny, image_x_lower, image_x_upper, image_y_lower, image_y_upper)
  require(grDevices)
  require(plotrix)
  color2D.matplot(res$stab_image[nrow(res$stab_image):1,], border=NA, xlab=NA, ylab=NA, axes=FALSE)
  #color2D.matplot(res$stab_image, border=NA, xlab=NA, ylab=NA, axes=FALSE)
  xlen <- length(res$xvec)
  ylen <- length(res$yvec)
  xarr <- c(1, as.integer((xlen+1)/4), as.integer((xlen+1)/2), as.integer((xlen+1)*3/4), xlen)
  yarr <- c(1, as.integer((ylen+1)/4), as.integer((ylen+1)/2), as.integer((ylen+1)*3/4), ylen)
  axis(1, at=xarr+0.5, labels=res$xvec[xarr])
  axis(2, at=yarr+0.5, labels=res$yvec[yarr])
}

plotStabilityImage1 <- function(centroids, theta, image_nx, image_ny, image_x_lower = 0, image_x_upper = 0, image_y_lower = 0, image_y_upper = 0){
  res <- stab_image <- make_stability_image(centroids, theta, image_nx, image_ny, image_x_lower, image_x_upper, image_y_lower, image_y_upper)
  require(graphics)
  image(res$xvec, res$yvec, res$stab_image[nrow(res$stab_image):1,])
}
