#'
#' correctBatchFSPY
#'
#' @description
#' Remove batch effect in FSPY object
#' @param object An FSPY object
#' @param batch vector. Batch covariate (only one batch allowed)
#' @param par.prior logical. TRUE indicates parametric adjustments will be used,
#'    FALSE indicates non-parametric adjustments will be used.
#' @param prior.plots logical. TRUE give prior plots with black as a kernel
#'    estimate of the empirical batch effect density and red as the parametric.
#' @param mean.only logical. FALSE If TRUE ComBat only corrects the mean of the batch
#'    effect (no scale adjustment)
#' @param ... Parameters passing to \code{\link[ComBat]{sva}} function
#'
#' @seealso \code{\link[BiocNeighbors]{findKNN}}
#'
#' @return An FSPY object after removing batch effect
#'
#' @importFrom sva ComBat
#' @export
#'
#'
correctBatchFSPY <- function(object, batch = NULL, par.prior = TRUE,
                             mean.only = TRUE, verbose = FALSE) {
  log.data <- object@log.data

  # correct batch effect using ComBat
  log.data.combat <- ComBat(dat = t(log.data), batch = batch,
                            par.prior = par.prior,
                            mean.only = mean.only, ...)

  object@log.data <- t(log.data.combat)
  return(object)
}

