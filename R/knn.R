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
runKNN <- function(object, knn = 30, BNPARAM = KmknnParam(),
                   knn.replace = F,
                   iter.max = 10,
                   knn.cluster = T,
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


#
# root.cells
#
addTrunk <- function(object) {
  trunk.network <- list()
  trunk.network$trunk.marker <- aggregate(object@log.data, list(cluster = object@meta.data$trunk.id), mean)
  rownames(trunk.network$trunk.marker) <- trunk.network$trunk.marker$cluster
  trunk.network$trunk.marker <- trunk.network$trunk.marker[, -1]
  trunk.network$trunk.dist <- stats::dist(trunk.network$trunk.marker, method = "euclidean")
  trunk.network$trunk.graph <- igraph::graph.adjacency(as.matrix(trunk.network$trunk.dist),
                                         mode = "undirected",
                                         weighted = TRUE)
  trunk.network$trunk.spanning.tree <- igraph::minimum.spanning.tree(trunk.network$trunk.graph)
  object@trunk.network <- trunk.network
  return(object)
}


#
# root.cells
#
addBranch <- function(object) {
  trunk.id <- table(object@meta.data$trunk.id)
  object@meta.data$branch.id <- 0
  branch.sub.graph <- vector("list", length(names(trunk.id)))
  names(branch.sub.graph) <- names(trunk.id)
  for (trunk.id.sub in names(trunk.id)) {
    cells.sub.idx <- which(object@meta.data$trunk.id == trunk.id.sub)
    sub.i <- igraph::graph.adjacency(object@network$knn.G[cells.sub.idx, cells.sub.idx], mode="undirected")
    sub.i <- simplify(sub.i)
    branch.sub.graph[[trunk.id.sub]] <- sub.i
    km.sub <- igraph::cluster_walktrap(sub.i)
    object@meta.data$branch.id[cells.sub.idx] <- paste0(trunk.id.sub, "-" , km.sub$membership)
  }
  branch.network <- list(branch.sub.graph = branch.sub.graph)
  branch.network$branch.marker <- aggregate(object@log.data, list(cluster = object@meta.data$branch.id), mean)
  rownames(branch.network$branch.marker) <- branch.network$branch.marker$cluster
  branch.network$branch.marker <- branch.network$branch.marker[, -1]
  branch.network$branch.dist <- stats::dist(branch.network$branch.marker, method = "euclidean")
  branch.network$branch.graph <- igraph::graph.adjacency(as.matrix(branch.network$branch.dist),
                                          mode = "undirected",
                                          weighted = TRUE)
  branch.network$branch.spanning.tree <- igraph::minimum.spanning.tree(branch.network$branch.graph)
  object@branch.network <- branch.network
  return(object)
}

#'
#' mclust
#'
#' @import mclust
#'
addMclust <- function(object) {
  mc.id <- Mclust(object@log.data)$classification

  object@meta.data$mc.id <- mc.id

  return(object)
}


slingshot <- function(object) {
  sds <- slingshot(object@log.data, object@meta.data$mc.id, start.clus = '1')

  plot(rd, col = cl, asp = 1)
  lines(sds, lwd = 3)

  return(object)
}


