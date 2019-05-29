#'
#' Calculate k-nearest neighbors of FSPY
#'
#' @name runKNN
#'
#' @description Calculates and stores a k-nearest neighbor graph based on Euclidean
#'    distance with (KMKNN) algorithm using log-transformed signaling matrix of
#'    flow cytometry data. The base function are base on \code{\link[BiocNeighbors]{findKNN}}.
#'
#' @param object an FSPY object
#' @param knn numeric. Number of k-nearest neighbors.
#' @param knn.replace logic. Whether to replace knn in FSPY object
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[BiocNeighbors]{findKNN}} function
#'
#' @seealso \code{\link[BiocNeighbors]{findKNN}}
#'
#' @return An FSPY object with knn, knn.index and knn.distance information.
#'
#' @import BiocNeighbors
#' @export
#'
#'
runKNN <- function(object,
                   knn = 30,
                   knn.replace = T,
                   verbose = F, ...) {

  if (isTRUE(object@knn > 0) & !(knn.replace)) {
    if (verbose) message(Sys.time(), " [INFO] Using knn in FSPY object: ", object@knn )
  } else if ( isTRUE(object@knn > 0) & (knn.replace) ) {
    if (verbose) message(Sys.time(), " [INFO] Using knn provided in this function: ", knn )
    object@knn <- knn
  } else {
    object@knn <- knn
  }

  if (verbose) message(paste0(Sys.time(), " [INFO] Calculating KNN " ) )
  fout <- findKNN(object@log.data, k = object@knn, ...)

  rownames(fout$index) <- object@meta.data$cell
  rownames(fout$distance) <- object@meta.data$cell

  object@knn = knn
  object@knn.index = fout$index
  object@knn.distance = fout$distance

  if (verbose) message(Sys.time(), " [INFO] Calculating KNN completed. ")
  return(object)
}

