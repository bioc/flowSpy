#'
#' Apply gating on FSPY object
#'
#' @param object S4 object. FSPY object
#' @param gating data.frame. gating paramete, the first column in gating is the marker name,
#'    the second column is the start paramter of gating cutoff, and the third column
#'    is the end paramter of gating cutoff.
#' @param initialization logical. Whether to initialize gating process
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
gatingFSPY <- function(object, gating = NULL, initialization = T, verbose = T) {
  if (is.null(gating)) {
    message(paste0(Sys.time(), " [WARNING] gating parameter is missing"))
    return(object)
  } else {
    if (!is.data.frame(gating)) {
      message(paste0(Sys.time(), " [WARNING] gating parameter must be a data frame"))
      return(object)
    } else {
      if (verbose) message( paste0(Sys.time(), " [INFO] start gating") )
      gate.data <- object@log.data
      if (initialization) {
        object@gating <- gating
      } else {
        object@gating <- rbind(object@gating, gating)
      }
      for (i in 1:nrow(object@gating)) {
        idx <- which(colnames(gate.data) == as.character(gating$marker[i]))
        gate.data <- gate.data[ which((gate.data[,idx] > as.numeric(gating$gs[i])) & (gate.data[,idx] < as.numeric(gating$ge[i])) ), ]
      }
      if (nrow(gate.data) <= 0) message(paste0(Sys.time(), " [WARNING] No cell remains after gating"))
      object@gate.data <- gate.data
      if (verbose) message( paste0(Sys.time(), " [INFO] Filteration of raw.meta.data") )
      object@raw.meta.data$gate = 0
      idx = match(rownames(gate.data), object@raw.meta.data$cell)
      object@raw.meta.data$gate[idx] = 1
      object@meta.data <- object@raw.meta.data[idx, ]
      return(object)
    }
  }

}




