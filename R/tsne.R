#'
#' Calculate t-Distributed Stochastic Neighbor Embedding in FSPY
#'
#' @param object an FSPY object
#' @param dims integer, Output dimensionality (default: 2)
#' @param initial_dims integer. the number of dimensions that should
#'    be retained in the initial PCA step (default: 50). See \code{\link[Rtsne]{Rtsne}}
#' @param perplexity numeric. Perplexity parameter. See \code{\link[Rtsne]{Rtsne}}
#' @param theta numeric. Speed/accuracy trade-off (increase for less accuracy),
#'    set to 0.0 for exact TSNE (default: 0.5). See \code{\link[Rtsne]{Rtsne}}
#' @param check_duplicates logical. Checks whether duplicates are present.
#'    It is best to make sure there are no duplicates present and set this
#'    option to FALSE, especially for large datasets (default: TRUE).
#'    See \code{\link[Rtsne]{Rtsne}}
#' @param verbose logical. Whether to print calculation progress.
#' @param pca,max_iter,is_distance,Y_init,pca_center,pca_scale See \code{\link[Rtsne]{Rtsne}}
#'
#' @import Rtsne
#'
#' @references
#'    Maaten, L. Van Der, 2014. Accelerating t-SNE using Tree-Based
#'    Algorithms. Journal of Machine Learning Research, 15, p.3221-3245.
#'
#'    van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing High-Dimensional
#'    Data Using t-SNE. Journal of Machine Learning Research, 9, pp.2579-2605.
#'
#' @export
#'
#' @examples
#'
#'
runTSNE <- function(object, dims = 2, initial_dims = 50, perplexity = 30,
                    theta = 0.5, check_duplicates = TRUE, pca = TRUE, max_iter = 1000,
                    verbose = TRUE, is_distance = FALSE, Y_init = NULL,
                    pca_center = TRUE, pca_scale = FALSE) {

  # tSNE calculation
  if (verbose) message(Sys.time(), " [INFO] Calculating tSNE.")
  tsne.obj <- Rtsne(as.matrix(object@log.data),
                    dims = dims, initial_dims = initial_dims, perplexity = perplexity,
                    theta = theta, check_duplicates = check_duplicates, pca = pca, max_iter = max_iter,
                    verbose = FALSE, is_distance = is_distance, Y_init = Y_init,
                    pca_center = pca_center, pca_scale = pca_scale)

  object@tsne.value <- tsne.obj$Y
  colnames(object@tsne.value) <- paste0("tSNE", 1:ncol(tsne.obj$Y))
  rownames(object@tsne.value) <- rownames(object@log.data)

  if (verbose) message(Sys.time(), " [INFO] Calculating tSNE completed. ")

  return(object)
}







