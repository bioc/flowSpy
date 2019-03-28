#'
#' Calculate principal components in FSPY
#'
#' @param object an FSPY object
#' @param center logical, a logical value indicating whether the variables
#'    should be shifted to be zero centered. Alternately, a vector
#'    of length equal the number of columns of x can be supplied.
#'    The value is passed to scale. See \code{\link[gmodels]{fast.prcomp}}
#' @param scale. logical, a logical value indicating whether the
#'    variables should be scaled to have unit variance before the
#'    analysis takes place. The default is FALSE for consistency
#'    with S, but in general scaling is advisable. Alternatively,
#'    a vector of length equal the number of columns of x can be supplied.
#'    The value is passed to scale. See \code{\link[gmodels]{fast.prcomp}}
#' @param verbose logical. Whether to print calculation progress.
#'
#' @importFrom gmodels fast.prcomp
#'
#' @examples
#'
#' @export
#'
runFastPCA <- function(object, center = FALSE, scale. = TRUE,  verbose = T) {
  # PCA calculation
  if (verbose) message(Sys.time(), " [INFO] Calculating PCA.")
  pca.obj <- fast.prcomp( t(object@log.data), retx = TRUE, center = center, scale. = scale.)

  object@pca.sdev <- pca.obj$sdev
  object@pca.load <- pca.obj$rotation
  object@pca.scores <- pca.obj$x

  if (verbose) message(Sys.time(), " [INFO] Calculating PCA completed. ")

  return(object)
}










