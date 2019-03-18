#' FSPY show method
#'
#' Prevents R from crashing by trying to print all slots of an FSPY object
#' if a returned object is not stored in a variable.
#'
#' @param object An FSPY object
#'
#' @exportMethod show
#'
#' @name show
#'
#' @docType methods
#'
#' @keywords internal
#'
#'
setMethod(
  f = "show",
  signature = "FSPY",
  definition = function(object) {
    cat(
      "FSPY Information:\n",
      "FSPY object:", nrow(object@raw.data), " cells \n",
      "Gating:", nrow(object@gating), " parameters \n",
      "FSPY object:", nrow(object@gate.data), " cells after gating \n"
    )
    invisible(NULL)
  }
)
