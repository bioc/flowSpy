#'
#' Update plot metadata of FSPY
#'
#' @name updatePlotMeta
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
#'
updatePlotMeta <- function(object, verbose = T) {
  plot.meta <- object@meta.data
  if (dim(object@log.data)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@log.data)
  }
  if (dim(object@pca.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@pca.value)
  }
  if (dim(object@tsne.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@tsne.value)
  }
  if (dim(object@dm@eigenvectors)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@dm@eigenvectors)
  }
  if (dim(object@umap.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@umap.value)
  }

  plot.meta <- as.data.frame(plot.meta)
  object@plot.meta <- plot.meta

  if (verbose) message(Sys.time(), " [INFO] Columns can be used in plot2D and plot3D: ", paste(colnames(plot.meta), collapse = " "))

  return(object)
}


#'
#' Fetching plot metadata of FSPY
#'
#' @name fetchPlotMeta
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
#'
fetchPlotMeta <- function(object, verbose = F) {

  object <- updatePlotMeta(object, verbose = verbose)

  return(object@plot.meta)
}



