#'
#' definition of root cells
#'
#' @name defRootCells
#'
#' @description definition of root cells
#'
#' @param object an FSPY object
#' @param root.cells vector. Cell name of the root cells
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @export
#'
#'
defRootCells <- function(object, root.cells = NULL, verbose = F) {
  if (length(object@root.cells) != 0) message(Sys.time(), " [INFO] root.cells in FSPY object exist, they will be replaced.")

  if (!is.vector(root.cells)) stop(Sys.time(), " [ERROR] root.cells must be a vector")

  if (is.character(root.cells)) {
    root.cells <- root.cells[root.cells %in% object@meta.data$cell]
  } else if (is.numeric(root.cells)) {
    root.cells <- object@meta.data$cell[object@meta.data$cluster.id %in% root.cells]
  } else {
    stop(Sys.time(), " [ERROR] invalid root.cells .")
  }

  ds.cells <- object@meta.data$cell[which(object@meta.data$dowsample == 1)]
  root.cells <- root.cells[root.cells %in% ds.cells]

  object@meta.data$is.root.cells <- 0
  object@meta.data$is.root.cells[match(root.cells, object@meta.data$cell)] <- 1
  if ( length(root.cells) == 0 ) {
    stop(Sys.time(), " [ERROR] root.cells are not in meta.data")
  } else {
    object@root.cells <- root.cells
  }

  if (verbose) message(Sys.time(), " [INFO] ", length(root.cells),  " cells will be added to root.cells .")

  return(object)
}

#'
#' definition of leaf cells
#'
#' @name defLeafCells
#' @description definition of root cells
#'
#' @param object an FSPY object
#' @param leaf.cells character or numeric. Cell name of the root cells or cluster.id of root.cells
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @export
#'
#'
defLeafCells <- function(object, leaf.cells = NULL, pseudotime.cutoff = 0, verbose = F) {
  if (length(object@leaf.cells) != 0) message(Sys.time(), " [INFO] leaf.cells in FSPY object exist, they will be replaced.")

  if (!is.vector(leaf.cells)) stop(Sys.time(), " [ERROR] leaf.cells must be a vector")

  if (is.character(leaf.cells)) {
    leaf.cells <- leaf.cells[leaf.cells %in% object@meta.data$cell]
  } else if (is.numeric(leaf.cells)) {
    leaf.cells <- object@meta.data$cell[object@meta.data$cluster.id %in% leaf.cells]
  } else {
    stop(Sys.time(), " [ERROR] invalid leaf.cells .")
  }

  ds.cells <- object@meta.data$cell[which(object@meta.data$dowsample == 1)]
  leaf.cells <- leaf.cells[leaf.cells %in% ds.cells]

  if (pseudotime.cutoff > 0) {
    if ( !all("pseudotime" %in% colnames(object@meta.data)) ) {
      warning(Sys.time(), " [WARNING] pseudotime is not in meta.data of FSPY, please run Pseudotime first.")
      pseudotime.cutoff = 0
    }
  }

  leaf.time <- object@meta.data$pseudotime[match(leaf.cells, object@meta.data$cell)]
  leaf.cells <- leaf.cells[which(leaf.time > pseudotime.cutoff )]

  object@meta.data$is.leaf.cells <- 0
  object@meta.data$is.leaf.cells[match(leaf.cells, object@meta.data$cell)] <- 1
  if ( length(leaf.cells) == 0 ) {
    stop(Sys.time(), " [ERROR] leaf.cells are not in meta.data")
  } else {
    object@leaf.cells <- leaf.cells
  }

  if (verbose) message(Sys.time(), " [INFO] ", length(leaf.cells),  " cells will be added to leaf.cells .")

  return(object)
}



#'
#' Calculation of Pseudotime
#'
#' @name runPseudotime
#'
#' @description calculation of Pseudotime based on KNN
#'
#' @param object An FSPY object
#' @param mode character. Specifies how igraph should interpret the supplied matrix.
#'    Possible values are: directed, undirected, upper, lower, max, min, plus.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to calculation function.
#'
#' @importFrom igraph graph.adjacency simplify distances
#' @return An FSPY object
#'
#' @export
#'
runPseudotime <- function(object, mode = "undirected", verbose = F, ...) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing.")

  if (dim(object@knn.index)[1] == 0) stop(Sys.time(), " [ERROR] KNN information is missing in FSPY, please run runKNN first.")

  if (verbose) message(Sys.time(), " [INFO] Calculating Pseudotime.")

  if ("pseudotime" %in% colnames(object@meta.data)) message(Sys.time(), " [INFO] Pseudotime exists in meta.data, it will be replaced.")

  knn.index <- object@knn.index

  adj <- matrix(0, nrow(knn.index), nrow(knn.index))
  rownames(adj) <- colnames(adj) <- rownames(knn.index)
  for(i in seq_len(ncol(knn.index))) {
    adj[i, rownames(knn.index)[knn.index[i,]]] <- 1
  }
  g <- igraph::graph.adjacency(adj, mode = mode, ... )
  # remove self loops
  g <- simplify(g)


  dist.all.path <- distances(g, v = as.character(object@root.cells))
  dist.all.path[which(is.infinite(dist.all.path))] <- max(dist.all.path[which(!is.infinite(dist.all.path))]) + 1
  pst <- colMeans(dist.all.path, na.rm = TRUE)
  pst <- ( pst - min(pst) )/ max( pst - min(pst) )

  object@meta.data$pseudotime <- 0
  object@meta.data$pseudotime[which(object@meta.data$dowsample == 1)] <- pst
  object@meta.data$traj.value <- 0
  object@meta.data$traj.value.log <- 0

  if (verbose) message(Sys.time(), " [INFO] Calculating Pseudotime completed.")

  return(object)
}









