#'
#' runCompensation
#'
#' @name runCompensation
#'
#' @param x An object of class flowFrame or flowSet.
#' @param spillover matrix, The spillover or compensation matrix.
#' @param ...  options to pass on to the compensate functions. See \code{\link[compensate]{flowCore}}
#'
#' @export
#'
#' @example
#'
#' # DO NOT RUN
#' if (F) {
#'
#' }
#'
runCompensation <- function(x, spillover = NULL, ...) {

  if (missing(x)) stop(Sys.time(), " [ERROR] input object is missing")

  if (is.null(spillover)) {
    if (is.null(data.1@description$SPILL)) {
      stop(Sys.time(), " [ERROR] compensation matrix is missing")
    } else {
      spillover <- x@description$SPILL
    }
  } else {
    if (!is.matrix(spillover)) {
      stop(Sys.time(), " [ERROR] spillover must be a matrix")
    }
  }

  x <- compensate(x, comp.mat)

  return(x)
}





