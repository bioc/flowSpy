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
#' @param knn.replace logic. Whether to replace knn in FSPY object
#' @param knn.cluster logic. Whether to run knn cluster.
#' @param verbose logic. Whether to print calculation progress.
#'
#'
#' @return An FSPY object
#'
#' @import BiocNeighbors
#' @importFrom igraph graph.adjacency cluster_walktrap
#' @export
#'
#' @examples
#'
#'
runKNN <- function(object, knn = 30, BNPARAM = KmknnParam(),
                   knn.replace = F,
                   knn.cluster = T,
                   verbose = T) {

  if (isTRUE(object@knn > 0) & !(knn.replace)) {
    if (verbose) message(paste0(Sys.time(), " [INFO] using knn in FSPY object: ", object@knn ) )
  } else if ( isTRUE(object@knn > 0) & (knn.replace) ) {
    if (verbose) message(paste0(Sys.time(), " [INFO] using knn provided in this function: ", knn ) )
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

  if (knn.cluster) {
    if (verbose) message(paste0(Sys.time(), " [INFO] start running knn.cluster. It will take some minutes if the dataset is too large. " ) )
    mat <- t(object@log.data)
    adj <- matrix(0, ncol(mat), ncol(mat))
    rownames(adj) <- colnames(adj) <- colnames(mat)
    for(i in seq_len(ncol(mat))) {
      adj[i,colnames(mat)[object@knn.index[i,]]] <- 1
    }
    g <- igraph::graph.adjacency(adj, mode="undirected")
    g <- simplify(g)
    ## identify communities
    km <- igraph::cluster_walktrap(g)
    object@meta.data$cluster.id <- km$membership
  } else {
    if (!"cluster.id" %in% colnames(object@meta.data)) {
      object@meta.data$cluster.id <- 0
    }
  }
  if (verbose) message(paste0(Sys.time(), " [INFO] Calculating KNN completed. "))
  return(object)
}



