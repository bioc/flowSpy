#'
#' definition of root cells
#'
#' @param object an FSPY object
#' @param root.cells vector, cell name of the root cells
#'
#' @export
#'
#' @examples
#'
#'
#'
defRootCells <- function(object, root.cells = NULL) {
  if (length(object@root.cells) != 0) warning(Sys.time(), " [WARNING] root.cells in FSPY object exist, they will be replaced.")

  if (!is.vector(root.cells)) stop(Sys.time(), " [ERROR] root.cell must be a vector")

  root.cells <- root.cells[root.cells %in% object@meta.data$cell]
  if ( length(root.cells) == 0 ) {
    stop(Sys.time(), " [ERROR] root.cell not in meta.data")
  } else {
    object@root.cells <- root.cells
  }

  return(object)
}


#'
#' Pseudotime
#'
Pseudotime <- function(object) {

  mat <- t(object@log.data)
  adj <- matrix(0, ncol(mat), ncol(mat))
  rownames(adj) <- colnames(adj) <- colnames(mat)
  for(i in seq_len(ncol(mat))) {
    adj[i, colnames(mat)[object@knn.index[i,]]] <- 1
  }
  g <- igraph::graph.adjacency(adj, mode="undirected")
  # remove self loops
  g <- simplify(g)


  dist.all.path <- distances(g, v = object@root.cells)
  pst <- colMeans(dist.all.path)

  object@meta.data$pseudotime <- pst

  return(pst)
}









