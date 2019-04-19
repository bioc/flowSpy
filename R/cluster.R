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
    if (verbose) message(Sys.time(), " [INFO] using knn in FSPY object: ", object@knn )
  } else if ( isTRUE(object@knn > 0) & (knn.replace) ) {
    if (verbose) message(Sys.time(), " [INFO] using knn provided in this function: ", knn )
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



