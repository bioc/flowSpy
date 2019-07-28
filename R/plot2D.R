#'
#' Visualization of 2D data of FSPY
#'
#' @name plot2D
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#' @param order.by vector. Order of color theme.
#' @param size numeric. Size of the dot
#' @param alpha numberic. Transparency (0-1) of the dot, default is 1.
#' @param category character. numeric or categorical
#' @param show.cluser.id logical. Whether to show cluster id in the plot.
#' @param show.cluser.id.size numeric. Size of the cluster id.
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
plot2D <- function(object,
                   item.use = c("PC_1", "PC_2"),
                   color.by = "stage",
                   order.by = NULL,
                   size = 1,
                   alpha = 1,
                   category = "categorical",
                   show.cluser.id = F,
                   show.cluser.id.size = 4,
                   main = "2D plot of FSPY",
                   plot.theme = theme_bw()) {

  # update and fetch plot meta information
  plot.meta <- fetchPlotMeta(object, verbose = F)
  idx <- match(c(color.by, item.use), colnames(object@log.data))
  idx <- idx[which(!is.na(idx))]
  if (length(idx) > 0) {
    sub <- as.data.frame(object@log.data[which(object@meta.data$dowsample == 1), idx])
    colnames(sub) <- colnames(object@log.data)[idx]
    plot.meta <- data.frame(plot.meta, sub)
  }

  # check item.use parameter in plot.meta data.frame
  if ( !all(item.use %in% colnames(plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  # check color.by parameter in plot.meta data.frame
  if ( !all(color.by %in% colnames(plot.meta)) ) stop(Sys.time(), " [ERROR] color.by is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if (length(item.use) < 2) stop(Sys.time(), " [ERROR] item.use is less than two elements.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " [WARNING] item.use has more than two elements. Only the first two will be used")
    item.use <- item.use[1:2]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " [WARNING] color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }

  item.use.idx <- match(item.use, colnames(plot.meta))
  color.by.idx <- match(color.by, colnames(plot.meta))

  plot.data <- data.frame(plot.x = plot.meta[, item.use.idx[1]],
                          plot.y = plot.meta[, item.use.idx[2]],
                          color.by = plot.meta[, color.by.idx])

  if ((length( unique(plot.data$color.by) ) > 256) & (category != "numeric")) {
    warning(Sys.time(), " [WARNING] color.by is categorical and has more than 256 elements. It will be used as numeric instead.")
    category = "numeric"
  }

  if (is.null(category)) {
    if (is.numeric(plot.data$color.by)) category="numeric" else category="categorical"
  }
  if (category == "categorical") {
    if (is.null(order.by)) {
      plot.data$color.by <- factor(plot.data$color.by)
    } else {
      plot.data$color.by <- factor(as.character(plot.data$color.by), levels = order.by)
    }
  } else if (category == "numeric") {
    if (!is.numeric(plot.data$color.by)) plot.data$color.by <- as.numeric(factor(plot.data$color.by))
  } else {
    warning(Sys.time(), " [WARNING] Unidentified parameters of category")
  }

  # plot
  gg <- ggplot(plot.data) + geom_point(aes(x=plot.x, y=plot.y, color = color.by), size = size, alpha = alpha)
  gg <- gg + plot.theme
  gg <- gg + labs(x = item.use[1], y = item.use[2], title = paste0(main))

  if (show.cluser.id & (category == "categorical")) {
     pos <- aggregate(  plot.data[, 1:2], list( pos = plot.data$color.by ), mean)

     for ( i in 1:length(pos$pos)) {
       gg <- gg + annotate(geom="text", x = pos$plot.x[i], y = pos$plot.y[i], label = pos$pos[i],
                           size = show.cluser.id.size)
     }
  }

  return(gg)

}

#'
#' Visualization violin plot of FSPY
#'
#' @name plotViolin
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#' @param order.by vector. Order of color theme.
#' @param size numeric. Size of the dot
#' @param alpha numberic. Transparency (0-1) of the dot, default is 1.
#' @param category character. numeric or categorical
#' @param show.cluser.id logical. Whether to show cluster id in the plot.
#' @param show.cluser.id.size numeric. Size of the cluster id.
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
plotViolin <- function(object,
                       marker,
                       color.by = "cluster.id",
                       order.by = NULL,
                       size = 1,
                       text.angle = 0,
                       main = "Violin plot FSPY",
                       plot.theme = theme_bw()) {

  # update plot meta information
  plot.meta <- fetchPlotMeta(object, verbose = F)

  if (missing(marker)) stop(Sys.time(), " [ERROR] marker is missing.")
  # check item.use parameter in plot.meta data.frame
  if (length(marker) > 1) {
    warning(Sys.time(), " [WARNING] marker has more than two elements. Only the first two will be used")
    marker <- marker[1]
  }
  if ( !all(marker %in% colnames(object@log.data)) ) stop(Sys.time(), " [ERROR] marker name is not correct")
  plot.meta <- data.frame(plot.meta, marker = object@log.data[which(fspy@meta.data$dowsample == 1), marker])

  # check color.by parameter in plot.meta data.frame
  if ( !all(color.by %in% colnames(plot.meta)) ) stop(Sys.time(), " [ERROR] color.by is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if (length(color.by) > 1) {
    warning(Sys.time(), " [WARNING] color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }
  color.by.idx <- match(color.by, colnames(plot.meta))

  plot.data <- data.frame(marker.by = plot.meta$marker,
                          color.by = plot.meta[, color.by.idx])

  if (length( unique(plot.data$color.by) ) > 256) {
    warning(Sys.time(), " [WARNING] color.by is categorical and has more than 256 elements. It will be used as numeric instead.")
  }

  if (is.null(order.by)) {
    plot.data$color.by <- factor(plot.data$color.by)
  } else {
    plot.data$color.by <- factor(as.character(plot.data$color.by), levels = order.by)
  }

  # plot
  gg <- ggplot(plot.data, aes(x = color.by, y= marker.by, fill = color.by)) + geom_violin(scale = "width")
  gg <- gg + plot.theme
  gg <- gg + stat_summary(fun.y=mean, geom="point", size = size, color="black")
  gg <- gg + labs(y = marker, x = color.by, title = paste0(main))
  gg <- gg + theme(axis.text.x = element_text(angle = text.angle, hjust = 1, vjust = 1))

  return(gg)

}




#'
#' Visualization pie plot of cluster data of FSPY
#'
#' @name plotPieCluster
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#' @param order.by vector. Order of color theme.
#' @param size numeric. Size of the dot
#' @param alpha numberic. Transparency (0-1) of the dot, default is 1.
#' @param category character. numeric or categorical
#' @param show.cluser.id logical. Whether to show cluster id in the plot.
#' @param show.cluser.id.size numeric. Size of the cluster id.
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
plotPieCluster <- function(object,
                           item.use = c("PC_1", "PC_2"),
                           cex.size = 1,
                           size.by.cell.number = T,
                           main = "2D pie plot of FSPY",
                           plot.theme = theme_bw()) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  if (is.null(object@network)) stop(Sys.time(), " [ERROR] network is missing, please run runCluster first!")
  if (length(unique(object@meta.data$stage)) <= 1) stop(Sys.time(), " [ERROR] plotPieTree only fits elements in stage over 2!")

  # update plot meta information
  plot.data <- fetchClustMeta(object, verbose = F)

  # check item.use parameter in cluster data.frame
  if ( !all(item.use %in% colnames(object@cluster)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if (length(item.use) < 2) stop(Sys.time(), " [ERROR] item.use is less than two elements.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " [WARNING] item.use has more than two elements. Only the first two will be used")
    item.use <- item.use[1:2]
  }
  item.use.idx <- match(item.use, colnames(object@cluster))

  plot.cols <- paste0(unique(object@meta.data$stage), ".percent")

  plot.data <- data.frame(plot.data,
                          pos.x = object@cluster[, item.use.idx[1]],
                          pos.y = object@cluster[, item.use.idx[2]])

  gg <- ggplot()
  if (size.by.cell.number) {
    gg <- gg + geom_scatterpie(aes(x = pos.x, y = pos.y, group = cluster, r = cell.number.percent*cex.size),
                               data = plot.data, cols = plot.cols, color=NA) + coord_equal()
  } else {
    gg <- gg + geom_scatterpie(aes(x = pos.x, y = pos.y, group = cluster, r = 0.1*cex.size),
                               data = plot.data, cols = plot.cols, color=NA) + coord_equal()
  }

  gg <- gg + plot.theme
  gg <- gg + labs(x = "", y = "", title = main)

  return(gg)

}

#'
#' Visualization of cluster data of FSPY
#'
#' @name plotCluster
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#' @param order.by vector. Order of color theme.
#' @param size numeric. Size of the dot
#' @param alpha numberic. Transparency (0-1) of the dot, default is 1.
#' @param category character. numeric or categorical
#' @param show.cluser.id logical. Whether to show cluster id in the plot.
#' @param show.cluser.id.size numeric. Size of the cluster id.
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
plotCluster <- function(object,
                        item.use = c("PC_1", "PC_2"),
                        color.by = "cluster",
                        size.by = "cell.number.percent",
                        order.by = NULL,
                        size = 1,
                        alpha = 1,
                        category = "categorical",
                        show.cluser.id = F,
                        show.cluser.id.size = 4,
                        main = "2D plot of cluster in FSPY",
                        plot.theme = theme_bw()) {

  # update plot meta information
  plot.meta.data <- fetchClustMeta(object, verbose = F)
  plot.meta.data <- data.frame(plot.meta.data, object@cluster)

  # check item.use parameter in plot.meta data.frame
  if ( !all(item.use %in% colnames(plot.meta.data)) ) stop(Sys.time(), " [ERROR] item.use is not in cluster data of FSPY, please run processingCluster first.")

  # check color.by parameter in plot.meta data.frame
  if ( !all(color.by %in% colnames(plot.meta.data)) ) stop(Sys.time(), " [ERROR] color.by is not in cluster data of FSPY, please run processingCluster first.")

  # check size.by parameter in plot.meta data.frame
  if ( !all(size.by %in% colnames(plot.meta.data)) ) stop(Sys.time(), " [ERROR] size.by is not in cluster data of FSPY, please run processingCluster first.")

  if (length(item.use) < 2) stop(Sys.time(), " [ERROR] item.use is less than two elements.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " [WARNING] item.use has more than two elements. Only the first two will be used")
    item.use <- item.use[1:2]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " [WARNING] color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }
  if (length(size.by) > 1) {
    warning(Sys.time(), " [WARNING] size.by has more than one elements. Only the first one will be used")
    size.by <- size.by[1]
  }

  item.use.idx <- match(item.use, colnames(plot.meta.data))
  color.by.idx <- match(color.by, colnames(plot.meta.data))
  size.by.idx <- match(size.by, colnames(plot.meta.data))


  plot.data <- data.frame(plot.x = plot.meta.data[, item.use.idx[1]],
                          plot.y = plot.meta.data[, item.use.idx[2]],
                          color.by = plot.meta.data[, color.by.idx],
                          size.by = plot.meta.data[, size.by.idx])

  if ((length( unique(plot.data$color.by) ) > 256) & (category != "numeric")) {
    warning(Sys.time(), " [WARNING] color.by is categorical and has more than 256 elements. It will be used as numeric instead.")
    category = "numeric"
  }

  if (is.null(category)) {
    if (is.numeric(plot.data$color.by)) category="numeric" else category="categorical"
  }
  if (category == "categorical") {
    if (is.null(order.by)) {
      plot.data$color.by <- factor(plot.data$color.by)
    } else {
      plot.data$color.by <- factor(as.character(plot.data$color.by), levels = order.by)
    }
  } else if (category == "numeric") {
    if (!is.numeric(plot.data$color.by)) plot.data$color.by <- as.numeric(factor(plot.data$color.by))
  } else {
    warning(Sys.time(), " [WARNING] Unidentified parameters of category")
  }

  # plot
  gg <- ggplot(plot.data) + geom_point(aes(x=plot.x, y=plot.y, color = color.by, size = size*size.by), alpha = alpha)
  gg <- gg + plot.theme
  gg <- gg + labs(x = item.use[1], y = item.use[2], title = paste0(main))

  if (show.cluser.id) {
    for ( i in 1:nrow(plot.data)) {
      gg <- gg + annotate(geom="text", x = plot.data$plot.x[i], y = plot.data$plot.y[i],
                          label = rownames(plot.data)[i],
                          size = show.cluser.id.size)
    }
  }

  return(gg)

}


#'
#' Visualization heatmap of cluster data of FSPY
#'
#' @name plotClusterHeatmap
#'
#' @param object An FSPY object
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#'
#' @import pheatmap
#'
#' @export
#'
#' @examples
#'
plotClusterHeatmap <- function(object,
                               color = colorRampPalette(c("blue","white","red"))(100),
                               scale = "row", ...) {

  # update plot meta information
  plot.meta.data <- fetchClustMeta(object, verbose = F)

  mat <- plot.meta.data[, object@markers]
  rownames(mat) <- plot.meta.data$cluster
  gg <- pheatmap(t(mat), color = color, scale = scale, border_color = NA, ...)

  return(gg)

}


#'
#' Visualization heatmap of data of FSPY
#'
#' @name plotHeatmap
#'
#' @param object An FSPY object
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#'
#' @import pheatmap
#'
#' @export
#'
#' @examples
#'
plotHeatmap <- function(object,
                        color = colorRampPalette(c("blue","white","red"))(100),
                        scale = "row",
                        downsize = 1000,
                        cluster_rows = F,
                        cluster_cols = F,
                        ...) {

  # update plot meta information
  plot.meta.data <- fetchPlotMeta(object, verbose = F)

  if (downsize > dim(plot.meta.data)[1]) {
    warning(Sys.time(), " [ERROR] Too large sample size of downsampling")
    downsize = dim(plot.meta.data)[1]
  }
  plot.meta.data <- plot.meta.data[sample(1:dim(plot.meta.data)[1], downsize), ]

  if (max(plot.meta.data$pseudotime) > 0) plot.meta.data <- plot.meta.data[order(plot.meta.data$pseudotime), ]

  mat <- object@log.data[match(plot.meta.data$cell, rownames(object@log.data)), ]
  gg <- pheatmap(t(mat), color = color, scale = scale, cluster_rows = cluster_rows,
                 cluster_cols = cluster_cols, border_color = NA, fontsize_col = 0.01, ...)

  return(gg)

}



















