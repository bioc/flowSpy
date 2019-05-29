#'
#' correctBatchFSPY
#'
#' @description
#' Remove batch effect in FSPY object
#' @param object An FSPY object
#' @param batch vector.
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
correctBatchFSPY <- function(object, batch = NULL, verbose = F) {


  return(object)
}

