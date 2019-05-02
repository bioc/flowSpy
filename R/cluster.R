#'
#' Specific Clustering Method Toolkits
#'
#' @name runCluster
#'
#' @description Compute a specific clustering using the combined flow
#'    cytometry data. "som" \code{\link[SOM]{flowSOM}}, "hclust" \code{\link[hclust]{stats}},
#'    "mclust" \code{\link[mclust]{mlcust}}, "kmeans" \code{\link[kmeans]{stats}} are
#'    provided.
#'
#' @param object an FSPY object
#' @param cluster.method character.
#' @param verbose logic. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @export
#'
#' @examples
#'
#'
#'
#'
runCluster <- function(object, cluster.method = "som", verbose = T, ...) {

  if (missing(object)) {
    stop(Sys.time(), " [ERROR] FSPY object is missing ")
  }

  if (cluster.method == "som") {
    object <- runSOM(object, verbose = verbose, ...)
  } else if (cluster.method == "hclust") {
    object <- runHclust(object, verbose = verbose, ...)
  } else if (cluster.method == "mclust") {
    object <- runMclust(object, verbose = verbose, ...)
  } else if (cluster.method == "kmeans") {
    object <- runKmeans(object, verbose = verbose, ...)
  } else {
    warning(Sys.time(), " [WARNING] Invalid cluster.method parameter ")
  }

  return(object)

}





#'
#' Calculate k-nearest neighbors of FSPY
#'
#' @name runKNN
#'
#' @description Calculates and stores a k-nearest neighbor graph based on Euclidean
#'    distance with (KMKNN) algorithm using log-transformed signaling matrix of
#'    flow cytometry data. The base function referes to \code{\link[BiocNeighbors]{findKNN}}.
#'
#' @param object an FSPY object
#' @param knn numeric, number of k-nearest neighbors
#' @param BNPARAM A BiocNeighborParam object, or NULL if BININDEX is supplied.
#'    See \code{\link[BiocNeighbors]{findKNN}}.
#' @param iter.max numeric. The maximum number of iterations allowed.
#' @param knn.replace logic. Whether to replace knn in FSPY object
#' @param knn.cluster logic. Whether to run knn cluster.
#' @param verbose logic. Whether to print calculation progress.
#'
#'
#' @return An FSPY object
#'
#' @import BiocNeighbors
#' @importFrom igraph graph.adjacency cluster_walktrap minimum.spanning.tree
#' @export
#'
#' @examples
#'
#'
runKNN <- function(object, cluster.method = c("som", "kmeans", "mclust", "hclust"),
                   list.params = list(),
                   verbose = T) {

  if (isTRUE(object@knn > 0) & !(knn.replace)) {
    if (verbose) message(Sys.time(), " [INFO] Using knn in FSPY object: ", object@knn )
  } else if ( isTRUE(object@knn > 0) & (knn.replace) ) {
    if (verbose) message(Sys.time(), " [INFO] Using knn provided in this function: ", knn )
    object@knn <- knn
  } else {
    object@knn <- knn
  }

  if (verbose) message(paste0(Sys.time(), " [INFO] Calculating KNN " ) )
  fout <- findKNN(object@log.data, k = object@knn, BNPARAM = BNPARAM)

  rownames(fout$index) <- object@meta.data$cell
  rownames(fout$distance) <- object@meta.data$cell

  object@knn = knn
  object@knn.index = fout$index
  object@knn.distance = fout$distance

  if (verbose) message(Sys.time(), " [INFO] Calculating KNN completed. ")
  return(object)
}


#'
#' runHclust
#'
#' @name runHclust
#'
#' @description Hierarchical cluster analysis on a set of dissimilarities
#'    and methods for analyzing it.
#'
#' @param object an FSPY object
#' @param hclust.method character
#' @param dist.method character
#' @param k numeric
#' @param verbose verbose
#'
#'
#' @importFrom stats hclust
#'
#' @export
#'
#' @examples
#'
runHclust <- function(object, k = 25,
                      hclust.method = "complete", dist.method = "euclidean",
                      verbose = T, ...) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Hclust.")

  d <- dist(object@log.data, method = dist.method)
  hc <- hclust(d, method = hclust.method)
  hc.tree <- cutree(hc, k = k)

  object@meta.data$hclust.id <- hc.tree

  object@hclust <- list(k = k,
                        d = d,
                        hc = hc,
                        value = aggregate(object@log.data, list(cluster = object@meta.data$kmeans.id), mean))

  if (verbose) message(Sys.time(), " [INFO] Calculating Hclust completed.")
  return(object)
}



#'
#' runKmeans
#'
#' @name runKmeans
#'
#' @description Perform k-means clustering on a data matrix.
#'
#' @param object  an FSPY object
#'
#' @return an FSPY object
#'
#' @export
#'
#' @examples
#'
#'
#'
runKmeans <- function(object, k = 25, iter.max = 10, nstart = 1,
                      algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                      trace=FALSE, scale = T, verbose = T) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Kmeans.")

  if (scale) kmeans.data <- scale(object@log.data) else kmeans.data = object@log.data

  kmeans.info <- kmeans(kmeans.data, centers = k, iter.max = iter.max, nstart = nstart,
                        algorithm = algorithm, trace=FALSE)

  object@meta.data$kmeans.id <- kmeans.info$cluster

  object@kmeans <- list(k = k,
                        centers = kmeans.info$centers)

  if (verbose) message(Sys.time(), " [INFO] Calculating Kmeans completed.")
  return(object)
}


#'
#' runMclust
#'
#' @name runMclust
#'
#' @description Model-based clustering based on parameterized finite Gaussian mixture models.
#'
#' @param object  an FSPY object
#' @return an FSPY object
#'
#' @export
#'
#' @import mclust
#'
#'
runMclust <- function(object) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Mclust.")

  mod <- Mclust(object@log.data)

  object@meta.data$mclust.id <- mod$classification

  object@mclust <- list(data = mod$data,
                        BIC = mod$BIC,
                        z =  mod$z)

  if (verbose) message(Sys.time(), " [INFO] Calculating Mclust completed.")
  return(object)
}



#'
#' calculation SOM in FSPY object
#'
#' @description Build a self-organizing map
#'
#' @param object  an FSPY object
#' @param xdim  Width of the grid.
#' @param ydim  Hight of the grid.
#' @param rlen  Number of times to loop over the training data for each MST
#' @param mst   Number of times to build an MST
#' @param alpha Start and end learning rate
#' @param radius Start and end radius
#' @param init  Initialize cluster centers in a non-random way
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev,
#'              4=cosine)
#' @param codes Cluster centers to start with
#' @param importance array with numeric values. Parameters will be scaled
#'                   according to importance
#' @param method the distance measure to be used. This must be one of "euclidean",
#'      "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'      Any unambiguous substring can be given. See \code{\link[stats]{dist}}
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return an FSPY object
#' @seealso \code{\link{BuildSOM}}
#'
#' @references This code is strongly based on the \code{\link[FlowSOM]{SOM}} function.
#'             Which is developed by Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2018).
#'
#' @import FlowSOM
#' @importFrom igraph graph.adjacency minimum.spanning.tree layout.kamada.kawai
#'
#' @export
#'
#' @examples
#'
runSOM <- function(object, xdim = 5, ydim = 5, rlen = 8, mst = 1,
                   alpha = c(0.05,  0.01), radius = 1, init = FALSE,
                   distf = 2, codes = NULL, importance = NULL,
                   method = "euclidean", verbose= T, ...) {

  if (verbose) message(Sys.time(), " [INFO] Calculating FlowSOM.")
  # flowSOM
  flowset <- as.matrix(object@log.data)
  flowsom <- FlowSOM::SOM(flowset,
                          xdim = xdim, ydim = ydim, rlen = rlen, mst = mst,
                          alpha = alpha[1], radius = radius,
                          init = init,
                          distf = distf, silent = verbose,
                          codes = codes, importance = importance)


  # generating som network
  object@meta.data$som.id <- flowsom$mapping[, 1]
  object@meta.data$som.value <- flowsom$mapping[, 2]
  object@som <- flowsom

  #object@som.network <- buildSOMnet(flowsom, object, method = method)

  if (verbose) message(Sys.time(), " [INFO] Calculating FlowSOM completed.")
  return(object)
}













