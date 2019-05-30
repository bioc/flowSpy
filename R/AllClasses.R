#' @import Matrix
#' @import ggplot2
#' @import gmodels
#' @import Rtsne
#' @import destiny
#' @import FlowSOM
#' @import BiocNeighbors
#' @import matrixStats
#' @import flowUtils
#' @import umap
#' @import prettydoc
NULL

#' flowSpy class
#'
#' @name FSPYclass
#'
#' @description  All information stored in FSPY object. You can use \code{creatFSPY} to
#'    create this object. In this package, most of the functions will use
#'    FSPY object as input, and return a modified FSPY obejct as well.
#'
#' @slot raw.data matrix. Raw signal data captured using flow cytometry.
#' @slot log.data matrix. Log-transfromed dataset of raw.data.
#' @slot meta.data data.frame. Meta data information, and colnames of "stage" and "cell" are required.
#' @slot markers vector. Markers used in the calculation of PCA, tSNE, destiny and umap.
#' @slot markers.idx vector. Index of markers used in the calculation of PCA, tSNE, destiny and umap.
#' @slot cell.name vector. Cell names of log data.
#' @slot knn numeric. Numbers of nearest neighbors
#' @slot knn.index,knn.distance matrix. Each row of the \code{knn.index} matrix corresponds to a point
#'     in \code{log.data} and contains the row indices in \code{log.data} that are its nearest neighbors.
#'     And each row of the \code{knn.distance} contains the distance of its nearest neighbors.
#' @slot som list. Store som network information calculated using \code{\link[FlowSOM]{FlowSOM}}.
#' @slot hclust list. Hierarchical cluster analysis on a merged flow cytometry data set. See \code{\link[hclust]{stats}}
#' @slot mclust list. Gaussian mixture modelling for classfication. See \code{\link[mclust]{mclust}}
#' @slot kmeans list. Perform k-means clustering on the merged flow cytometry data set. See \code{\link[kmeans]{stats}}
#' @slot pca.sdev,pca.value,pca.scores Pca information of FSPY object which are
#'     generated from \code{\link[gmodels]{fast.prcomp}}.
#' @slot tsne.value matrix. tSNE coordinates information. See \code{\link[Rtsne]{Rtsne}}.
#' @slot dm DiffusionMap object. Diffusion map calculated by \code{\link[destiny]{DiffusionMap}}
#' @slot umap.value matrix umap coordinates information calculated using \code{\link[umap]{umap}}.
#' @slot root.cells vector, Names of root cells, which can be modified by \code{defRootCells}.
#' @slot leaf.cells vector. Names of leaf cells, which can be modified by \code{defLeafCells}.
#' @slot network list. Network stored in the calculation of trajectory and pseudotime.
#' @slot walk list. Random forward and backward walk between \code{root.cells} and \code{leaf.cells}.
#' @slot diff.tree list. Differentiation tree of all cells.
#' @slot diff.traj list. Differentiation trajectory all cells.
#' @slot plot.meta data.frame. Plot meta information for \code{plot2D} or \code{plot3D}.
#' @slot add.meta list. Additional meta information of FSPY object.
#'
#' @importClassesFrom destiny DiffusionMap DPT
#'
#' @exportClass FSPY
#'
#' @aliases FSPYclass
#'
#'
FSPY <- methods::setClass("FSPY", slots = c(
  raw.data = "matrix",
  log.data = "matrix",
  meta.data = "data.frame",
  markers = "vector",
  markers.idx = "vector",
  cell.name = "vector",

  # KNN
  knn = "numeric",
  knn.index = "matrix",
  knn.distance = "matrix",

  # som network
  som = "list",

  # hclust network
  hclust = "list",

  # mclust network
  mclust = "list",

  # kmeans network
  kmeans = "list",

  # pca information
  pca.sdev = "vector",
  pca.value = "matrix",
  pca.scores = "matrix",

  # tsne information
  tsne.value = "matrix",

  # diffusion map information
  dm = c("DiffusionMap", NULL),

  # umap information
  umap.value = "matrix",

  # run for pseudotime
  root.cells = "vector",
  leaf.cells = "vector",
  network = "list",

  # trajectory analysis
  walk = "list",
  diff.tree = "list",
  diff.traj = "list",

  # for visualization
  plot.meta = "data.frame",
  add.meta = "list"
  )
)

#' create an FSPY object
#'
#' @description
#' This function include all information used in flow cytometry data
#' processing.
#'
#' @name createFSPY
#'
#' @param raw.data matrix. Raw data read from FCS file after perform preprocessing.
#' @param markers vector. Detailed marker information in the gate of flow cytometer.
#' @param meta.data data.frame. Raw metadata of each cell. Columns "cell" and "stage" are required.
#' @param log.transform logical. Whether to log transformed raw.data. If FALSE, it's better
#'    to perform transformation using \code{\link[transformation]{flowCore}} before creating FSPY
#'    object. flowSpy only provide log transforma parameter. If you need to using truncateTransform,
#'    scaleTransform, linearTransform, quadraticTransform and lnTransform, see \code{flowCore} for more
#'    information.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... paramters pass to \code[correctBatchFSPY] function.
#'
#' @export
#'
#' @return An FSPY object with raw.data and markers and meta.data
#'
#' @examples
#'
#'
#'
createFSPY <- function(raw.data, markers, meta.data,
                       batch = NULL, batch.correct = F,
                       log.transform = T, verbose = F, ...) {
  # QC of cells
  if (missing(raw.data)) stop(Sys.time(), " [ERROR] raw.data is required")
  if (!is.matrix(raw.data)) {
    warning(Sys.time(), " [WARNING] raw.data must be a matrix")
    raw.data <- as.matrix(raw.data)
  }
  if (verbose) message(Sys.time(), " [INFO] Number of cells in processing: ", dim(raw.data)[1])

  # QC of metadata
  if (missing(meta.data)) stop(Sys.time(), " [ERROR] meta.data must be a data.frame")
  if (!is.data.frame(meta.data)) {
    warning(Sys.time(), " [WARNING] meta.data must be a data.frame")
    meta.data <- as.matrix(meta.data)
  }

  if (!all(c("cell", "stage") %in% colnames(meta.data))) {
    stop(Sys.time(), " [ERROR] cell and stage information must be provided in meta.data")
  }

  if (nrow(raw.data) != nrow(meta.data)) {
    stop(Sys.time(), " [ERROR] cell number in raw.data is not equal to that in meta.data")
  } else {
    if (verbose) message(Sys.time(), " [INFO] rownames of meta.data and raw.data will be named using column cell")
    rownames(raw.data) = meta.data$cell
    rownames(meta.data) = meta.data$cell
  }

  # load index of markers of FCS
  if (missing(markers)) stop(Sys.time(), " [ERROR] markers is missing")
  if (!is.vector(markers)) {
    warning(Sys.time(), " [WARNING] markers must be a vector")
    markers <- as.vector(markers)
  }
  markers.idx <- match(markers, colnames(raw.data))
  if (verbose) message(Sys.time(), " [INFO] Index of markers in processing")
  if (any(is.na(markers.idx))) {
    sub.markers <- markers[which(is.na(markers.idx))]
    warning(Sys.time(), " [WARNING] ", sub.markers, " not existes in colnames of raw.data. It will be removed. ")

    markers <- markers[which(!is.na(markers.idx))]
    markers.idx <- markers.idx[which(!is.na(markers.idx))]
  }

  # Create an FSPY object
  if (verbose) message(Sys.time(), " [INFO] Creating FSPY object.")
  object <- new("FSPY", raw.data = raw.data, meta.data = meta.data,
                markers = markers, markers.idx = markers.idx)

  # log transfromed using base 10
  if (log.transform) {
    if (verbose) message(Sys.time(), " [INFO] Log transformed ")
    object@log.data <- log10(abs(raw.data[, markers.idx]) + 1)
  } else {
    if (verbose) message(Sys.time(), " [INFO] No log transformed ")
    object@log.data <- raw.data[, markers.idx]
  }

  # correcting batch effect
  if (batch.correct) {
    if (is.null(batch)) {
      warning(Sys.time(), " [WARNING] batch must be provided when batch.correct is TRUE ")
    } else {
      object <- correctBatchFSPY(object, batch = batch, ...)
    }
  }

  object@plot.meta <- data.frame(row.names = object@meta.data$cell)
  object@meta.data$is.root.cells <- 0
  object@meta.data$is.leaf.cells <- 0

  if (verbose) message(Sys.time(), " [INFO] Build FSPY object succeed ")
  return(object)
}













