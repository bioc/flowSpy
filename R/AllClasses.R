#' @import Matrix
#' @import ggplot2
#' @import ggthemes
#' @import gmodels
#' @import Rtsne
#' @import destiny
#' @import FlowSOM
#' @import matrixStats
#' @import flowUtils
NULL

#' flowSpy class
#'
#' All information stored in FSPY object, and you can use creatFSPY to
#' create this object. In this package, most of the function will use
#' FSPY object as input, and return a modified SPTY obejct as well.
#'
#' @slot raw.data matrix, raw signal data captured using flow cytometry.
#' @slot raw.meta.data data.frame, raw meta data information from the experiment.
#' @slot log.data matrix, data set log transfromed from the raw.data.
#' @slot gate.data matrix, data set used in the calculation of PCA, tSNE and destiny.
#' @slot meta.data data.frame, meta data information after gating.
#' @slot markers vector, markers used in the calculation of PCA, tSNE and destiny.
#' @slot markers.idx vector, index of markers used in the calculation of PCA, tSNE and destiny.
#' @slot gating data.frame, gating parameter stored in the FSPY object.
#' @slot pca.sdev,pca.load,pca.score pca information of FSPY object which generated from \code{\link[gmodels]{fast.prcomp}}.
#' @slot tsne.value data.frame, tSNE coordinates information of FSPY object.
#' @slot dm data.frame, destiny information of FSPY object.
#' @slot som list, som network of FSPY object.
#' @slot cluster data.frame, cluster in FSPY.
#' @slot root.cell vector, names of root cells.
#' @slot leaf.cell vector, names of leaf cells.
#' @slot pseudotime data.frame, pseudotime of all cells.
#' @slot walk list, random forward and backward walk between tips and cells.
#' @slot diff.tree list, differentiation tree of FSPY object.
#' @slot diff.traj list, differentiation trajectory of FSPY object.
#' @slot network list, network of FSPY object.
#' @slot sankey list, sankey diagrams of FSPY object.
#'
#' @importClassesFrom destiny DiffusionMap
#'
#' @exportClass FSPY
#'
#' @aliases FSPYclass
#'
#' @name FSPY
#'
FSPY <- setClass("FSPY", slots = c(
  raw.data = "matrix",
  raw.meta.data = "data.frame",
  log.data = "matrix",
  gate.data = "matrix",
  meta.data = "data.frame",
  markers = "vector",
  markers.idx = "vector",
  gating = "data.frame",
  # pca information
  pca.sdev = "vector",
  pca.load = "matrix",
  pca.score = "matrix",
  # tsne information
  tsne.value = "data.frame",
  # diffusion map information
  dm = c("DiffusionMap", NULL),
  # som information
  som = "list",

  # run for improved function
  cluster = "data.frame",
  root.cell = "vector",
  leaf.cell = "vector",
  pseudotime = "data.frame",
  walk = "list",
  diff.tree = "list",
  diff.traj = "list",
  network = "list",
  sankey = "list"
  )
)

#' create an FSPY object
#'
#' This function include all information used in flow cytometry data
#' processing.
#'
#' @param raw.data matrix. raw data read from FCS file.
#' @param markers vector. detailed marker information in the gate of flow cytometer
#' @param raw.meta.data data.frame. raw metadata of each cell. column "cell" and "stage" are required
#' @param gating data.frame. gating paramete, the first column in gating is the marker name,
#'    the second column is the start paramter of gating cutoff, and the third column
#'    is the end paramter of gating cutoff.
#' @param log.transformed logical. Whether to log transformed the raw data. If not, it's better
#'    to set the print parameter as logical.
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @examples
#'
#'
#' @export
#'
createFSPY <- function(raw.data, markers, raw.meta.data, gating = NULL, log.transformed = T, verbose = T) {
  # QC of cells
  if (missing(raw.data) | !is.matrix(raw.data)) stop(paste0(Sys.time(), " [ERROR] raw.data must be a matrix"))
  if (verbose) message(paste0(Sys.time(), " [INFO] Number of cells in processing: ", dim(raw.data)[1]))

  # QC of metadata
  if (missing(raw.meta.data) | !is.data.frame(raw.meta.data)) stop(paste0(Sys.time(), " [ERROR] raw.meta.data must be a data.frame"))
  if (nrow(raw.data) != nrow(raw.meta.data)) stop(paste0(Sys.time(), " [ERROR] cell number in raw.data is not equal to that in raw.meta.data")) else rownames(raw.data) = raw.meta.data$cell

  # load index of markers of FCS
  if (missing(markers) | !is.vector(markers)) stop(paste0(Sys.time(), " [ERROR] markers must be a vector"))
  markers.idx <- match(markers, colnames(raw.data))
  if (verbose) message(paste0(Sys.time(), " [INFO] Index of markers in processing"))

  # Create an FSPY object
  if (verbose) message(paste0(Sys.time(), " [INFO] Creating FSPY object."))
  object <- new("FSPY", raw.data = raw.data, raw.meta.data = raw.meta.data,
                markers = markers, markers.idx = markers.idx)

  # log transfromed using base 10
  if (log.transformed) {
    tag <- unlist(lapply(1:nrow(raw.data), function(x) sum(raw.data[x, markers.idx] < 1)  ))
    if (sum(tag) > 0) message(paste0(Sys.time(), " [WARNING] Fluorescence intensity lower than 1 will be filted"))
    if (verbose) message(paste0(Sys.time(), " [INFO] Log transformed "))
    object@log.data <- log10(raw.data[which(tag == 0), markers.idx])
  } else {
    if (verbose) message(paste0(Sys.time(), " [INFO] No log transformed "))
    object@log.data <- raw.data[, markers.idx]
  }

  # gating
  if (!is.null(gating)) object@gating <- gating
  object <- gatingFSPY(object, gating = gating, verbose = verbose)

  if (verbose) message(paste0(Sys.time(), " [INFO] Build FSPY object succeed "))
  return(object)
}


























