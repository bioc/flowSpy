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
                   item.use = c("PC1", "PC2"),
                   color.by = "stage",
                   order.by = NULL,
                   size = 1,
                   alpha = 1,
                   category = "categorical",
                   show.cluser.id = F,
                   show.cluser.id.size = 4,
                   main = "2D plot of FSPY",
                   plot.theme = theme_bw()) {

  # update plot meta information
  object <- updatePlotMeta(object, verbose = F)

  # check item.use parameter in plot.meta data.frame
  if ( !all(item.use %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  # check color.by parameter in plot.meta data.frame
  if ( !all(color.by %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] color.by is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if (length(item.use) < 2) stop(Sys.time(), " [ERROR] item.use is less than two elements.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " [WARNING] item.use has more than two elements. Only the first two will be used")
    item.use <- item.use[1:2]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " [WARNING] color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }

  item.use.idx <- match(item.use, colnames(object@plot.meta))
  color.by.idx <- match(color.by, colnames(object@plot.meta))

  plot.data <- data.frame(plot.x = object@plot.meta[, item.use.idx[1]],
                          plot.y = object@plot.meta[, item.use.idx[2]],
                          color.by = object@plot.meta[, color.by.idx])

  if ((length( unique(plot.data$color.by) ) > 50) & (category != "numeric")) {
    warning(Sys.time(), " [WARNING] color.by is categorical and has more than 50 elements. It will be used as numeric instead.")
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
#' Visualization pie plot of 2D data of FSPY
#'
#' @name plotPie2D
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
plotPie2D <- function(object,
                      item.use = c("PC_1", "PC_2"),
                      cex.size = 1,
                      size.by.cell.number = T,
                      main = "2D pie plot of FSPY",
                      plot.theme = theme_bw()) {

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



