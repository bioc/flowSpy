#'
#' Update plot metadata of FSPY
#'
#' @name updatePlotMeta
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
#'
updatePlotMeta <- function(object, verbose = T) {
  plot.meta <- object@meta.data
  if (dim(object@log.data)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@log.data)
  }
  if (dim(object@pca.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@pca.value)
  }
  if (dim(object@tsne.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@tsne.value)
  }
  if (dim(object@dm@eigenvectors)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@dm@eigenvectors)
  }
  if (dim(object@umap.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@umap.value)
  }

  plot.meta <- as.data.frame(plot.meta)
  object@plot.meta <- plot.meta

  if (verbose) message(Sys.time(), " [INFO] Columns can be used in plot2D and plot3D: ", paste(colnames(plot.meta), collapse = " "))

  return(object)
}


#'
#' Fetching plot metadata of FSPY
#'
#' @name fetchPlotMeta
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
#'
fetchPlotMeta <- function(object, verbose = F) {

  object <- updatePlotMeta(object, verbose = verbose)

  return(object@plot.meta)
}

#'
#' Fetching cellls of FSPY
#'
#' @name fetchCell
#'
#' @param object An FSPY object
#' @param logical.connect character. "and" or "or"
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Paramters to pass to limitation
#'
#' @export
#'
#'
fetchCell <- function(object, logical.connect = "or", verbose = F, ... ) {

  object <- updatePlotMeta(object, verbose = verbose)

  param.list <- list(...)
  plot.meta <- object@plot.meta

  cell.left <- NULL

  if (length(param.list) > 0) {
    param.list <- param.list[names(param.list) %in% colnames(plot.meta)]
    for (i in 1:length(param.list)) {
      sub <- param.list[[i]]
      sub.name <- names(param.list)[i]
      if (sub.name %in% c("cell", "stage", "is.root.cells", "is.leaf.cells")) {
        cell.sub <- plot.meta[plot.meta[, sub.name] %in% sub, "cell"]
      } else if (grepl(".id$", sub.name)) {
        cell.sub <- plot.meta[plot.meta[, sub.name] %in% sub, "cell"]
      } else if (is.numeric(sub)) {
        if (length(sub) == 1) {
          cell.sub <- plot.meta[which(plot.meta[, sub.name] > sub[1]), "cell"]
        } else {
          cell.sub <- plot.meta[which((plot.meta[, sub.name] > sub[1]) & (plot.meta[, sub.name] < sub[2])), "cell"]
        }
      } else {
        cell.sub <- NULL
      }
      if (logical.connect == "and") {
        cell.left <- intersect(cell.left, cell.sub)
      } else if (logical.connect == "or") {
        cell.left <- union(cell.left, cell.sub)
      } else {
        stop(Sys.time(), " [ERROR] Unidentified logical.connect")
      }
    }
  }

  return(cell.left)
}




#'
#' constraintMatrix
#'
#' @description
#' constraint FCS data by a provid cutoff
#'
#' @param x matrix
#' @param cutoff numeric. Cutoff of the constraint value
#' @param markers character. Markers used in the calculation of constraint model.
#' @param method character. the distance measure to be used.
#'    This must be one of "euclidean", "maximum", "manhattan",
#'    "canberra", "binary" or "minkowski".
#'
#' @export
#'
constraintMatrix <- function(x, cutoff = 0.99, markers = NULL, method = "euclidean") {

  if (!is.numeric(x)) stop(Sys.time(), " [ERROR] x must be a matrix ")

  if (is.null(markers)) markers <- colnames(x)
  if (!all(markers %in% colnames(x))) stop(Sys.time(), " [ERROR] markers must belong to the colnames of x ")

  sub <- abs(x[, markers])
  if (length(markers) > 1) {
    sub.mean <- colMeans(sub)
    d <- sapply(1:nrow(sub), function(aa) dist(rbind(sub[aa,], sub.mean), method = method))
  } else {
    sub.mean <- mean(sub)
    d <- sapply(1:length(sub), function(aa) dist(rbind(sub[aa], sub.mean), method = method))
  }

  filter.mat <- x[order(d), ]
  filter.mat <- filter.mat[1:floor(cutoff*dim(filter.mat)[1]), ]

  return(filter.mat)
}








