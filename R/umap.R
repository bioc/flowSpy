#'
#' Calculating UMAP
#'
#' @name runUMAP
#'
#' @description
#' Calculate Uniform Manifold Approximation and Projection in FSPY
#'
#' @param object an FSPY object
#' @param umap.config object of class umap.config. See \code{\link[umap]{umap}}.
#' @param dim numeric. Dim of umap, you can also change it in umap.config.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Options to pass on to the \code{\link[umap]{umap}} function
#'
#' @import umap
#' @seealso \code{\link[umap]{umap}}
#'
#' @export
#'
runUMAP <- function(object, umap.config = umap.defaults, dim = 2, verbose = F, ...) {
  if (verbose) message(Sys.time(), " [INFO] Calculating Umap.")

  umap.config$n_neighbors <- object@knn
  umap.config$n_components <- dim
  umap.out <- umap(object@log.data, config = umap.config, ...)
  object@umap.value <- umap.out$layout
  colnames(object@umap.value) <- paste0("UMAP_", 1:ncol(umap.out$layout))
  rownames(object@umap.value) <- rownames(object@log.data)

  if (verbose) message(Sys.time(), " [INFO] Calculating Umap.")
  return(object)
}
