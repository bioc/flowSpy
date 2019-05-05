
#'
#' Walk between root cells and leaf cells
#'
#' @name runWalk
#'
#' @description Walk between root cells and leaf cells
#'
#' @param object An FSPY object
#' @param mode character. Specifies how igraph should interpret the supplied matrix.
#'    Possible values are: directed, undirected, upper, lower, max, min, plus.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to calculation function.
#'
#' @importFrom igraph graph.adjacency simplify shortest_paths
#' @return An FSPY object
#'
#' @export
#'
runWalk <- function(object, mode = "undirected",
                    max.run = 20, weighted = T,
                    verbose = T, ...) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing.")

  if (dim(object@knn.index)[1] == 0) stop(Sys.time(), " [ERROR] KNN information is missing in FSPY, please run runKNN first.")

  if (verbose) message(Sys.time(), " [INFO] Calculating walk between root.cells and leaf.cells .")

  if (!"pseudotime" %in% colnames(object@meta.data)) stop(Sys.time(), " [INFO] Pseudotime exists in meta.data, it will be replaced.")

  mat <- t(object@log.data)
  adj <- matrix(0, ncol(mat), ncol(mat))
  rownames(adj) <- colnames(adj) <- colnames(mat)
  pseudotime <- object@meta.data$pseudotime
  for(i in seq_len(ncol(mat))) {
    idx <- object@knn.index[i,][ pseudotime[object@knn.index[i,]] > pseudotime[i]  ]
    adj[i, colnames(mat)[idx]] <- 1
  }
  g <- igraph::graph.adjacency(adj, mode = mode, ...)
  # remove self loops
  g <- simplify(g)

  root.cells <- object@root.cells
  if (length(root.cells) > max.run ) {
    root.cells <- as.character(sample(root.cells, max.run))
  }
  walk <- lapply(root.cells, function(x) shortest_paths(g, from = x, to = object@leaf.cells)$vpath )

  cell.info <- object@meta.data$cell[unlist(walk)]
  cell.info <- as.data.frame(table(cell.info))

  object@meta.data$traj.value <- cell.info$Freq[match( object@meta.data$cell, cell.info$cell.info)] / max.run
  object@meta.data$traj.value[object@root.cells] = 0
  object@meta.data$traj.value[object@leaf.cells] = 0

  object@walk <- list(max.run = max.run,
                      walk = walk)

  object@meta.data$traj.value.log <- log10(object@meta.data$traj.value + 1)

  if (verbose) message(Sys.time(), " [INFO] Calculating walk completed.")

  return(object)
}





