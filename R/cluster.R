#'
#' Specific Clustering Method Toolkits
#'
#' @name runCluster
#'
#' @description Compute a specific clustering using the combined flow
#'    cytometry data. "som" \code{\link[flowSOM]{SOM}}, "hclust" \code{\link[stats]{hclust}},
#'    "mclust" \code{\link[mclust]{mlcust}}, "kmeans" \code{\link[stats]{kmeans}} are
#'    provided.
#'
#' @param object an FSPY object
#' @param cluster.method character. Four clustering method are provided: som, hclust, mclust, kmeans.
#' @param verbose logic. Whether to print calculation progress.
#' @param ... options to pass on to the clustering functions.
#'
#' @seealso Four clustering methods are provided: \code{\link[flowSOM]{SOM}}, \code{\link[stats]{hclust}},
#'    \code{\link[mclust]{mlcust}}, \code{\link[stats]{kmeans}}. You can use \code{runSOM}, \code{runHclust},
#'    \code{runMclust} and \code{runKmeans} to run clustering respectively.
#'
#' @return An FSPY object with cluster.id in meta.data
#'
#'
#'
runCluster <- function(object, cluster.method = "som", verbose = T, ...) {

  if (missing(object)) {
    stop(Sys.time(), " [ERROR] FSPY object is missing ")
  }

  if (cluster.method == "som") {
    object <- runSOM(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$som.id
  } else if (cluster.method == "hclust") {
    object <- runHclust(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$hclust.id
  } else if (cluster.method == "mclust") {
    object <- runMclust(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$mclust.id
  } else if (cluster.method == "kmeans") {
    object <- runKmeans(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$kmeans.id
  } else {
    warning(Sys.time(), " [WARNING] Invalid cluster.method parameter ")
  }

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
#' @param hclust.method character or a function. The agglomeration method to be used.
#'    This should be one of "ward.D", "ward.D2", "single", "complete", "average",
#'    "mcquitty", "median" or "centroid". Or you can specify an equation as input, for example
#'    \code{function(x) hclust(x,method = 'ward.D2')}.
#' @param dist.method character or a function. The distance measure to be used.
#'    This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary"
#'    or "minkowski". Or you can specify an equation as input, for example
#'    \code{function(x) as.dist((1-cor(t(x)))/2)}.
#' @param k numeric. The number of clusters.
#' @param verbose logical. Whether to print calculation progress.
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link[stats]{dist}}
#'
#' @importFrom stats hclust dist
#'
#' @return An FSPY object with hclust.id in FSPY object
#' @export
#'
#'
runHclust <- function(object, k = 25,
                      hclust.method = "complete", dist.method = "euclidean",
                      verbose = T) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Hclust.")

  # check dist parameters
  if (is.character(dist.method)) {
    d <- stats::dist(object@log.data, method = dist.method)
  } else if (is.function(dist.method)) {
    d <- dist.method(object@log.data)
  } else {
    warning(Sys.time(), " [WARNING] Invalid dist.method parameter.")
    d <- stats::dist(object@log.data)
  }

  # check hclust parameters
  if (is.character(hclust.method)) {
    hc <- stats::hclust(d, method = hclust.method)
  } else if (is.function()) {
    hc <- dist.method(d)
  } else {
    warning(Sys.time(), " [WARNING] Invalid hclust.method parameter.")
    hc <- stats::hclust(d)
  }


  hc.tree <- cutree(hc, k = k)

  object@meta.data$hclust.id <- hc.tree

  object@hclust <- list(k = k,
                        d = d,
                        hc = hc,
                        value = aggregate(object@log.data, list(cluster = object@meta.data$hclust.id), mean))

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
#' @param k numeric. The number of clusters.
#' @param iter.max numeric. The maximum number of iterations allowed.
#' @param nstart numeric. If k is a number, how many random sets should be chosen.
#' @param algorithm numeric.
#' @param trace logical or integer number.
#' @param scale logical. Whether to use scales data in hclust.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[stats]{kmeans}} function
#'
#' @return an FSPY object with kmeans.id in meta.data
#'
#' @seealso \code{\link[stats]{kmeans}}
#'
#' @importFrom stats kmeans
#' @export
#'
#'
#'
#'
runKmeans <- function(object, k = 25, iter.max = 10, nstart = 1,
                      algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                      trace=FALSE, scale = F, verbose = T, ...) {

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
#'    This function is based on \code{\link[mclust]{Mclust}}.
#'
#' @param object  an FSPY object
#' @param scale logical. Whether to use scaled data in Mclust.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[mclust]{Mclust}} function
#'
#'
#' @return an FSPY object with mclust.id in meta.data
#'
#' @seealso \code{\link[mclust]{Mclust}}
#'
#' @export
#'
#' @importFrom mclust Mclust
#'
#'
runMclust <- function(object, scale = F,
                      verbose = T, ...) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Mclust.")

  if (scale) mclust.data <- scale(object@log.data) else mclust.data = object@log.data

  mod <- Mclust(mclust.data, ...)

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
#' @param ... Parameters passing to \code{\link[FlowSOM]{SOM}} function
#'
#' @return an FSPY object with som.id in FSPY object
#' @seealso \code{\link{BuildSOM}}
#'
#' @references This code is strongly based on the \code{\link[FlowSOM]{SOM}} function.
#'             Which is developed by Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2018).
#'
#' @importFrom FlowSOM SOM
#'
#' @seealso \code{\link[FlowSOM]{SOM}}
#'
#' @export
#'
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
                          codes = codes, importance = importance,
                          ...)


  # generating som network
  object@meta.data$som.id <- flowsom$mapping[, 1]
  object@meta.data$som.value <- flowsom$mapping[, 2]
  object@som <- flowsom

  #object@som.network <- buildSOMnet(flowsom, object, method = method)

  if (verbose) message(Sys.time(), " [INFO] Calculating FlowSOM completed.")
  return(object)
}













