#'
#' plot Pseudotime density of FSPY
#'
#' @name plotPseudotimeDensity
#'
#' @param object an FSPY object
#' @param color.by character.
#' @param main character. Title of the plot
#' @param adjust numeric. A multiplicate bandwidth adjustment.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @export
#'
plotPseudotimeDensity <- function(object, color.by = "stage",
                                  main = "Density of pseudotime",
                                  adjust = 0.5,
                                  plot.theme = theme_bw()) {
  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  object <- updatePlotMeta(object, verbose = F)

  # checking items
  if ( !all("pseudotime" %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] pseudotime is not in plot.meta of FSPY, please run Pseudotime first.")

  if ( !all(color.by %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  item.use.idx <- match("pseudotime", colnames(object@plot.meta))
  color.by.idx <- match(color.by, colnames(object@plot.meta))

  plot.data <- data.frame(pseudotime = object@plot.meta[, item.use.idx[1]],
                          color.by = object@plot.meta[, color.by.idx])

  if (length(unique(plot.data$color.by)) > 50) {
    message(Sys.time(), " [INFO] color.by is a numeric vector and has over 50 elements")
    plot.data$color.by <- 1
  }
  if (!is.factor(plot.data$color.by)) plot.data$color.by <- as.factor(as.character(plot.data$color.by))

  gg <- ggplot(plot.data, aes(x=pseudotime, colour = color.by))
  gg <- gg + geom_density(adjust = adjust)
  gg <- gg + plot.theme

  gg <- gg + labs(title = paste0(main))

  return(gg)
}


#'
#' plotPseudotimeTraj
#'
#' @name plotPseudotimeTraj
#'
#' @param object An FSPY object
#' @param markers character. Markers used in the calculation progress
#' @param cutoff numeric. Cutoff of trajectory value
#' @param size numeric. Size of the dot
#' @param alpha numeric. Transparency (0-1) of the dot, default is 1.
#' @param print.curve logical. Whether to perform curve fitting
#' @param var.cols logical. Whether to plot stage
#' @param plot.theme themes from \code{ggplot2}
#'
#' @importFrom stats predict
#'
#' @export
#'
plotPseudotimeTraj <- function(object,
                               cutoff = -1,
                               markers = NULL,
                               size = 0.5,
                               alpha = 0.6,
                               print.curve = T,
                               var.cols = F,
                               plot.theme = theme_bw()) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  object <- updatePlotMeta(object, verbose = F)

  # checking items
  if ( !all("pseudotime" %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] pseudotime is not in plot.meta of FSPY, please run Pseudotime first.")

  # checking items
  if (cutoff > 0) {
    if ( !all(c("traj.value","traj.value.log") %in% colnames(object@plot.meta)) ) {
      message(Sys.time(), " [INFO] traj.value is not in plot.meta of FSPY, please run runWalk first.")
    }
  }
  if (is.null(markers)) markers <- object@markers

  plot.data <- NULL
  plot.meta <- object@plot.meta
  for (i in 1:length(markers)) {
    sub <- data.frame(Pseudotime = plot.meta$pseudotime,
                      IsRoot = plot.meta$is.root.cells,
                      IsLeaf = plot.meta$is.leaf.cells,
                      Marker = markers[i],
                      Signal = plot.meta[, colnames(plot.meta) %in% markers[i]],
                      Stage = plot.meta$stage,
                      TrajValue = plot.meta$traj.value,
                      LogTrajValue = plot.meta$traj.value.log)
    idx <- which( (sub$LogTrajValue > cutoff)  )
    if (sum(sub$IsRoot == 1) > 0) {
      idx.2 <- which(sub$IsRoot == 1)
    } else {
      idx.2 <- NULL
    }
    idx <- union(idx, idx.2)
    if (sum(sub$IsLeaf == 1) > 0) {
      idx.3 <- which(sub$IsLeaf == 1)
    } else {
      idx.3 <- NULL
    }
    idx <- union(idx, idx.3)
    sub <- sub[idx, ]
    plot.data <- rbind(plot.data, sub)
  }

  gg <- ggplot(plot.data, aes(x=Pseudotime, y=Signal, color = Pseudotime)) + geom_point(size = size, alpha = alpha)
  gg <- gg + plot.theme

  if (var.cols) {
    gg <- gg + facet_grid(rows = vars(Marker), cols = vars(Stage))
  } else {
    gg <- gg + facet_grid(rows = vars(Marker))
  }

  if (print.curve) {
    gg <- gg + geom_smooth(color="black", method="loess", se = F)
  }

  return(gg)
}


#'
#' plotMarkerDensity
#'
#' @name plotMarkerDensity
#'
#' @param object An FSPY object
#' @param markers character. Markers used in the calculation progress
#' @param cutoff numeric. Cutoff of trajectory value
#' @param adjust numeric. Transparency (0-1) of the dot, default is 1.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @importFrom stats predict
#'
#' @export
#'
plotMarkerDensity <- function(object,
                              cutoff = -1,
                              markers = NULL,
                              adjust = 0.5,
                              plot.theme = theme_bw()) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  object <- updatePlotMeta(object, verbose = F)

  # checking items
  if (cutoff > 0) {
    if ( !all(c("traj.value","traj.value.log") %in% colnames(object@plot.meta)) ) {
      message(Sys.time(), " [INFO] traj.value is not in plot.meta of FSPY, please run runWalk first.")
    }
  }
  if (is.null(markers)) markers <- object@markers

  plot.data <- NULL
  plot.meta <- object@plot.meta
  for (i in 1:length(markers)) {
    sub <- data.frame(Pseudotime = plot.meta$pseudotime,
                      IsRoot = plot.meta$is.root.cells,
                      IsLeaf = plot.meta$is.leaf.cells,
                      Marker = markers[i],
                      Signal = plot.meta[, markers[i]],
                      Stage = plot.meta$stage,
                      TrajValue = plot.meta$traj.value,
                      LogTrajValue = plot.meta$traj.value.log)
    idx <- which( (sub$LogTrajValue > cutoff)  )
    if (sum(sub$IsRoot == 1) > 0) {
      idx.2 <- which(sub$IsRoot == 1)
    } else {
      idx.2 <- NULL
    }
    idx <- union(idx, idx.2)
    if (sum(sub$IsLeaf == 1) > 0) {
      idx.3 <- which(sub$IsLeaf == 1)
    } else {
      idx.3 <- NULL
    }
    idx <- union(idx, idx.3)
    sub <- sub[idx, ]
    plot.data <- rbind(plot.data, sub)
  }

  gg <- ggplot(plot.data, aes(x=Signal, color = Marker)) + geom_density(adjust = adjust)
  gg <- gg + plot.theme
  gg <- gg + facet_grid(rows = vars(Stage), cols = vars(Marker))
  return(gg)
}


#'
#' plot MST of FSPY
#'
#' @name plotTree
#'
#' @param object an FSPY object
#' @param cex.size numeric. size cex of the dot
#' @param color.by numeric. size color theme of the dot
#' @param size.by numeric. size theme of the dot
#' @param show.node.name logical. whether to show node name
#' @param color.theme character. Library of color theme
#'
#' @export
#' @importFrom igraph layout_as_tree layout.kamada.kawai as_data_frame
#'
plotTree <- function(object,
                     cex.size = 1,
                     color.by = "cell.percent",
                     size.by = "cell.percent",
                     as.tree = F,
                     root.id = NULL,
                     show.node.name = F,
                     color.theme = NULL) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  if (is.null(object@network)) stop(Sys.time(), " [ERROR] network is missing, please run runCluster first!")

  mst <- object@network$mst

  # update plot meta information
  plot.data <- fetchPlotMeta(object, verbose = F)

  cell.count <- table(plot.data[, match(c("cluster.id", "stage"), colnames(plot.data)) ])
  cell.percent <- sapply(1:dim(cell.count)[2], function(x) cell.count[,x]/sum(cell.count[,x]) )
  colnames(cell.percent) <- paste0(colnames(cell.count), ".percent")
  cell.percent.stage <- t(sapply(1:dim(cell.count)[1], function(x) cell.count[x,]/sum(cell.count[x,]) ))
  colnames(cell.percent.stage) <- paste0(colnames(cell.count), ".percent.stage")

  for (i in length(plot.data):1) {
    if (!is.numeric(plot.data[, i])) {
      plot.data <- plot.data[, -i]
    }
  }

  node.attr <- aggregate(plot.data, list(id = plot.data$cluster.id), mean)
  node.attr <- cbind(node.attr, cell.percent, cell.percent.stage)
  node.attr$cell.percent <- as.numeric(rowSums(cell.count))

  edge.attr <- igraph::as_data_frame(mst)

  ##### layout
  if (as.tree) {
    if (is.null(root.id)) {
      root.id = node.attr$id[which(node.attr$pseudotime == min(node.attr$pseudotime))]
    }
    l <- igraph::layout_as_tree(mst, root = root.id)
  } else {
    l <- igraph::layout.kamada.kawai(mst)
  }
  colnames(l) <- c("pos.x", "pos.y")
  node.attr <- cbind(node.attr, l)

  size.by.idx <- match(size.by, colnames(node.attr))
  color.by.idx <- match(color.by, colnames(node.attr))

  edge.attr$from.x <- node.attr$pos.x[match(edge.attr$from, node.attr$id)]
  edge.attr$from.y <- node.attr$pos.y[match(edge.attr$from, node.attr$id)]
  edge.attr$to.x <- node.attr$pos.x[match(edge.attr$to, node.attr$id)]
  edge.attr$to.y <- node.attr$pos.y[match(edge.attr$to, node.attr$id)]

  color.tree <- node.attr[, color.by.idx]
  size.tree <- node.attr[, size.by.idx]

  gg <- ggplot()
  gg <- gg + geom_segment(mapping = aes(x = edge.attr$from.x, y = edge.attr$from.y, xend = edge.attr$to.x, yend = edge.attr$to.y))
  gg <- gg + geom_point(mapping = aes(x = node.attr$pos.x, y = node.attr$pos.y, color = color.tree, size = size.tree))
  gg <- gg + scale_size(range = c(0, 6) * cex.size)

  if (show.node.name) gg <- gg + geom_text(aes(x = node.attr$pos.x, y = node.attr$pos.y, label = node.attr$id ), check_overlap = TRUE, size = 3 * cex.size)
  gg <- gg + theme_void() + coord_flip() + scale_y_reverse()
  gg <- gg + labs(x = "", y = "", title = paste0("Tree plot, color.by: ", color.by, ", size.by: ", size.by))

  return(gg)

}




