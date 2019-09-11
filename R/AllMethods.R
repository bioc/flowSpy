#' FSPY show method
#'
#' Prevents R from crashing by trying to print all slots of an FSPY object
#' if a returned object is not stored in a variable.
#'
#' @param object An FSPY object
#'
#' @exportMethod show
#' @export
#'
#' @name show
#'
#' @docType methods
#'
#' @keywords internal
#'
#' @return Cell number of FSPY object
#'
#' @examples
#' if (FALSE) {
#'   fspy
#' }
#'
setMethod(
  f = "show",
  signature = "FSPY",
  definition = function(object) {
    cat(
      "FSPY Information:\n",
      "Input cell number:", nrow(object@raw.data), " cells \n",
      "Enroll marker number:", ncol(object@log.data), " markers \n",
      "Cells after downsampling:", sum(object@meta.data$dowsample), " markers \n"
    )
    invisible(NULL)
  }
)
