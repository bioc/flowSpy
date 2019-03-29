#'
#' Update plot metadata of FSPY
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
#'
updatePlotMeta <- function(object, verbose = T) {
  plot.meta <- NULL
  if (dim(object@meta.data)[1] > 0) {
    
  }

  plot.meta <- cbind(object@meta.data, object@log.data, object@pca.load,
                     object@tsne.value, object@dm@eigenvectors, object@umap.layout)
  plot.meta <- as.data.frame(plot.meta)
  object@plot.meta <- plot.meta

  if (verbose) message(Sys.time(), " [INFO] Columns can be used in plot2D and plot3D: ", paste(colnames(plot.meta), collapse = " "))

  return(object)
}






