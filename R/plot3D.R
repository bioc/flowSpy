#'
#' Visualization of 3D data of FSPY
#'
#' @name plot3D
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 3D plot, axes x and y and z must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#' @param size numeric. size of the dot
#' @param angle numberic. angle of the plot
#' @param scale.y numeric. scale of y axis related to x- and z axis
#' @param main character. title of the plot
#' @param category character. numeric or categorical
#' @param plot.theme themes from \code{\link[ggthemes]{ggthemes}}
#' @param color.theme character. Library of color theme
#'
#' @import scatterplot3d
#'
#' @export
#'
#' @examples
#'
plot3D <- function(object,
                   item.use = c("PC1", "PC2", "PC3"),
                   color.by = "stage",
                   order.by = NULL,
                   size = 1,
                   angle = 60,
                   scale.y = 0.8,
                   category = NULL,
                   main = "3D plot of FSPY",
                   plot.theme = theme_base(),
                   color.theme = NULL) {

  object <- updatePlotMeta(object, verbose = F)

  if ( !all(item.use %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if ( !all(color.by %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if (length(item.use) < 3) stop(Sys.time(), " [ERROR] item.use is less than two characters.")
  if (length(item.use) > 3) {
    warning(Sys.time(), " [WARNING] item.use is more than two characters. Only the first two will be used")
    item.use <- item.use[1:3]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " [WARNING] color.by is more than one characters. Only the first one will be used")
    color.by <- color.by[1]
  }

  item.use.idx <- match(item.use, colnames(object@plot.meta))
  color.by.idx <- match(color.by, colnames(object@plot.meta))

  plot.data <- data.frame(plot.x = object@plot.meta[, item.use.idx[1]],
                          plot.y = object@plot.meta[, item.use.idx[2]],
                          plot.z = object@plot.meta[, item.use.idx[3]],
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
    plot.data$color.by.3d <- factor(plot.data$color.by, labels = rainbow(length(levels(plot.data$color.by))))
  } else if (category == "numeric") {
    if (!is.numeric(plot.data$color.by)) plot.data$color.by <- as.numeric(factor(plot.data$color.by))
    color.lib <- rainbow(102)
    plot.data$color.by.sd <- plot.data$color.by - min(plot.data$color.by)
    plot.data$color.by.3d <- color.lib[ ceiling( plot.data$color.by.sd/max(plot.data$color.by.sd) * 100 ) + 1 ]
  } else {
    warning(Sys.time(), " [WARNING] Unidentified parameters of category")
  }


  scatterplot3d(x = plot.data$plot.x, y = plot.data$plot.y, z = plot.data$plot.z,
                color = plot.data$color.by.3d,
                pch = 16, cex.symbols = size,
                scale.y = scale.y, angle = angle,
                xlab = item.use[1], ylab = item.use[2], zlab = item.use[3],
                main = main,
                col.axis = "#444444", col.grid = "#CCCCCC")


}
