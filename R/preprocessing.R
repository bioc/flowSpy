#'
#' runCompensation
#'
#' @description Run compensation by applying a spillover matrix.
#'    See \code{\link[flowCore]{compensate}} for more information.
#'
#' @name runCompensation
#'
#' @param x An object of class \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}.
#' @param spillover matrix, The spillover or compensation matrix.
#' @param ...  options to pass on to the further arguments. See \code{\link[flowCore]{compensate}}
#'
#' @seealso \code{\link[flowCore]{compensate}}
#'
#' @importFrom flowCore compensate
#' @export
#'
#' @examples
#'
#' # DO NOT RUN
#' if (F) {
#'   fileName <- system.file("extdata","D2.fcs",
#'                           package="flowSpy")
#'
#'   # read FACs data
#'   fcs.data <- read.FCS(filename = fileName)
#'
#'   # run compensation
#'   fcs.data <- runCompensation(fcs.data)
#'
#' }
#'
runCompensation <- function(x, spillover = NULL, ...) {

  if (missing(x)) stop(Sys.time(), " [ERROR] input object is missing")

  if (is.null(spillover)) {
    if (is.null(x@description$SPILL)) {
      stop(Sys.time(), " [ERROR] compensation matrix in x is missing")
    } else {
      spillover <- x@description$SPILL
    }
  } else {
    if (!is.matrix(spillover)) {
      stop(Sys.time(), " [ERROR] spillover must be a matrix")
    }
  }

  x <- flowCore::compensate(x, spillover, ...)

  return(x)
}





