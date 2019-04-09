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
#' @param size numeric. size of the dot
#' @param alpha numberic. transparency (0-1) of the dot, default is 0.2.
#' @param category character. numeric or categorical
#' @param main character. title of the plot
#' @param plot.theme themes from \code{\link[ggthemes]{ggthemes}}
#' @param color.theme character. Library of color theme
#'
#' @import ggplot2
#' @import ggthemes
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
                   category = NULL,
                   main = "2D plot of FSPY",
                   plot.theme = theme_base(),
                   trajectory = "som",
                   color.theme = NULL) {

  object <- updatePlotMeta(object, verbose = F)

  if ( !all(item.use %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if ( !all(color.by %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if (length(item.use) < 2) stop(Sys.time(), " [ERROR] item.use is less than two characters.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " [WARNING] item.use is more than two characters. Only the first two will be used")
    item.use <- item.use[1:2]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " [WARNING] color.by is more than one characters. Only the first one will be used")
    color.by <- color.by[1]
  }

  item.use.idx <- match(item.use, colnames(object@plot.meta))
  color.by.idx <- match(color.by, colnames(object@plot.meta))

  plot.data <- data.frame(plot.x = object@plot.meta[, item.use.idx[1]],
                          plot.y = object@plot.meta[, item.use.idx[2]],
                          color.by = object@plot.meta[, color.by.idx])

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
  gg <- ggplot() + geom_point(aes(x=plot.data$plot.x, y=plot.data$plot.y, color = plot.data$color.by), size = size, alpha = alpha)

  if (!is.null(color.theme)) {
    if (category == "categorical") {
      gg <- gg + scale_colour_gradientn(colours = color.theme)
    } else if (category == "numeric") {
      gg <- gg + scale_color_manual(values = color.theme)
    }
  }

  gg <- gg + theme_base()
  gg <- gg + labs(x = item.use[1], y = item.use[2], title = paste0(main))



  if (trajectory == "som") {
    if ("som.id" %in% colnames(object@meta.data)) {
      plot.data$som.id <- object@meta.data$som.id
      traj.layout.node <- aggregate(plot.data[, 1:2], list(som.id = plot.data$som.id), mean)


      som.net <- object@som.network
      traj.layout.edge <- igraph::as_data_frame(som.net$graph, what="edges")
      traj.layout.edge$from.x <- traj.layout.node$plot.x[match(traj.layout.edge$from, traj.layout.node$som.id)]
      traj.layout.edge$from.y <- traj.layout.node$plot.y[match(traj.layout.edge$from, traj.layout.node$som.id)]
      traj.layout.edge$to.x <- traj.layout.node$plot.x[match(traj.layout.edge$to, traj.layout.node$som.id)]
      traj.layout.edge$to.y <- traj.layout.node$plot.y[match(traj.layout.edge$to, traj.layout.node$som.id)]

      gg <- gg + geom_segment(aes(x = traj.layout.edge$from.x, y = traj.layout.edge$from.y, xend = traj.layout.edge$to.x, yend = traj.layout.edge$to.y))
      gg <- gg + geom_point(mapping = aes(x = traj.layout.node$plot.x, y = traj.layout.node$plot.y))

    }
  }

  return(gg)

}



