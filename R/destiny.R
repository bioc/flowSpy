#'
#' Calculate diffusion map in FSPY
#'
#' @param object an FSPY object
#' @param sigma.use numeric. Diffusion scale parameter of the Gaussian kernel. One of '\code{local}',
#'     '\code{global}', a \code{\link[base]{numeric}} global sigma or a Sigmas object.
#'     When choosing '\code{global}', a global sigma will be calculated using find_sigmas
#'     (See \code{\link[destiny]{find_sigmas}}). A larger sigma might be necessary if the eigenvalues can not
#'    be found because of a singularity in the matrix. See \code{\link[destiny]{DiffusionMap}}
#' @param distance Distance measurement method applied to data or a distance matrix/dist.
#'    For the allowed values, see \code{\link[destiny]{find_knn}}
#' @param density.norm logical. If TRUE, use density normalisation. See \code{\link[destiny]{DiffusionMap}}
#' @param verbose logical. Whether to print calculation progress.
#' @param ... options to pass on to the DiffusionMap function
#'
#' @seealso
#'
#' @import destiny
#'
#' @export
#'
runDiffusionMap <- function(object, sigma.use = NULL,
                            distance=c("euclidean", "cosine", "rankcor"),
                            density.norm = TRUE,  verbose=T,
                            ...) {

  dm.data <- as.matrix(object@log.data)

  if (verbose) message(paste0(Sys.time(), " [INFO] Calculating Diffusion Map."))
  # Figure out sigma
  # this function refered to URD calcDM function.
  if (is.null(sigma.use)) {
    sigma.use <- find_sigmas(dm.data, verbose=F)@optimal_sigma
    if (verbose) message(paste0(Sys.time(), " [INFO] Destiny determined an optimal global sigma: ", round(sigma.use, digits=3)))
  } else if (is.numeric(sigma.use)) {
    if (verbose) message(paste0(Sys.time(), " [INFO] Using provided global sigma: ", round(sigma.use, digits=3)))
  } else if (sigma.use == "local") {
    if (verbose) message(paste0(Sys.time(), " [INFO] Using local sigma "))
  } else {
    sigma.use <- find_sigmas(dm.data, verbose=F)@optimal_sigma
    warning(paste0(Sys.time(), " [WARNING] Invalid sigma value. Using an optimal global sigma instead."))
  }
  # Calculate the Diffusion Map
  dm.obj <- DiffusionMap(dm.data, sigma=sigma.use, k=object@knn, density_norm = density.norm, distance=distance[1], ...)

  rownames(dm.obj@eigenvectors) <- rownames(dm.data)
  rownames(dm.obj@transitions) <- rownames(dm.data)
  colnames(dm.obj@transitions) <- rownames(dm.data)

  # Load diffusion map into the Dropseq object
  object@dm <- dm.obj

  if (verbose) message(paste0(Sys.time(), " [INFO] Calculating Diffusion Map completed"))

  return(object)
}

