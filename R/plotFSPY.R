
#'
#' plot with gate
#'
#' @param object S4 object. FSPY object
#' @param plot.markers vector included 2 charaters. For example c("CD34", "CD43").
#' @param plot.type character. Plot style for flow cytometry data. default: dot.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of raw.meta.data, or it can be just "density" (the default value).
#' @param size numeric. size of the dot
#' @param alpha numberic. transparency (0-1) of the dot, default is 0.2.
#' @param color.theme character. Library of color theme
#' @param gated logical. Whether to show only gated data.
#' @param show.gate logical. Whether to show gated cells. This parameter will only take function when
#'     \code{color.by} is "density"
#'
#' @examples
#'
#' @export
#'
plotGATE <- function(object, plot.markers, plot.type = c("dot", "mesh", "dotmesh", "dotellipse","ellipse"),
                     color.by = "density",
                     size = 1,
                     alpha = 0.2,
                     color.theme = c("#666666", "#f2de00", "#8bd129", "#33CCFF", "#0066CC", "#CC66FF", "#FF3300"),
                     gated = T, show.gate = F) {
  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  if (missing(plot.markers)) stop(Sys.time(), " [ERROR] plot.markers are missing")
  if (!is.vector(plot.markers)) stop(Sys.time(), " [ERROR] plot.markers must be a vector, e.g. c(\"CD34\", \"CD43\")")
  if (length(plot.markers) != 2) stop(Sys.time(), " [ERROR] plot.markers must be a vector, e.g. c(\"CD34\", \"CD43\")")

  if (gated) {
    gg.data <- object@gate.data
    gg.data.plot <- data.frame(object@meta.data,
                               plot.x = gg.data[, which(colnames(gg.data) == plot.markers[1])],
                               plot.y = gg.data[, which(colnames(gg.data) == plot.markers[2])])

  } else {
    gg.data <- object@log.data
    gg.data.plot <- data.frame(object@raw.meta.data,
                               plot.x = gg.data[, which(colnames(gg.data) == plot.markers[1])],
                               plot.y = gg.data[, which(colnames(gg.data) == plot.markers[2])])

  }
  if (color.by == "density") {
    if (show.gate) gg.data.plot$color.by <- as.factor(gg.data.plot$gate) else gg.data.plot$color.by <- factor(0)
  } else {
    idx <- which(colnames(gg.data.plot) == color.by)
    if (is.factor(gg.data.plot[, idx])) {
      gg.data.plot$color.by <- gg.data.plot[, idx]
    } else {
      gg.data.plot$color.by <- factor(gg.data.plot[, idx])
    }

  }

  gg <- ggplot(gg.data.plot, aes(x=plot.x, y=plot.y, colour = color.by))
  if (plot.type[1] == "dot") {
    gg <- gg + geom_point(alpha = alpha, size = size )
  } else if (plot.type[1] == "mesh") {
    gg <- gg + geom_density_2d()
  } else if (plot.type[1] == "dotmesh") {
    gg <- gg + geom_point(alpha = alpha, size = size ) + geom_density_2d(aes(colour = color.by))
  } else if (plot.type[1] == "dotellipse") {
    gg <- gg + geom_point(alpha = alpha, size = size ) + stat_ellipse()
  } else if (plot.type[1] == "ellipse") {
    gg <- gg + stat_ellipse()
  }

  gg <- gg + scale_color_manual(values=color.theme)
  gg <- gg + theme_base()
  gg <- gg + labs(x = plot.markers[1], y = plot.markers[2], title = paste0(plot.type, " plot"))

  return(gg)
}


#'
#' Visualization of 2D data
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of raw.meta.data, or it can be just "density" (the default value).
#' @param size numeric. size of the dot
#' @param alpha numberic. transparency (0-1) of the dot, default is 0.2.
#' @param main character. title of the plot
#' @param color.theme character. Library of color theme
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
                   size = 1,
                   alpha = 1,
                   main = "2D plot of FSPY",
                   color.theme = c("#f2de00", "#8bd129", "#33CCFF", "#0066CC", "#CC66FF", "#FF3300")) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")

  # checking items
  colname.info <- c(colnames(object@pca.load),
                    colnames(object@tsne.value),
                    colnames(object@gate.data),
                    colnames(object@meta.data),
                    colnames(object@dm@eigenvectors),
                    colnames(object@pseudotime))

  if ( (sum(!item.use %in% colname.info) > 0) | (length(item.use) != 2) ) stop(paste0(Sys.time(), " [ERROR] invalided colname.info paramter "))

  plot.data <- data.frame(object@meta.data)

  # check paramter for plot.x
  if (item.use[1] %in% colnames(object@pca.load) ) {
    idx <- which(colnames(object@pca.load) == item.use[1])
    plot.data$plot.x <- object@pca.load[, idx]
  } else if (item.use[1] %in% colnames(object@tsne.value) ) {
    idx <- which(colnames(object@tsne.value) == item.use[1])
    plot.data$plot.x <- object@tsne.value[, idx]
  } else if (item.use[1] %in% colnames(object@gate.data) ) {
    idx <- which(colnames(object@gate.data) == item.use[1])
    plot.data$plot.x <- object@gate.data[, idx]
  } else if (item.use[1] %in% colnames(object@meta.data) ) {
    idx <- which(colnames(object@meta.data) == item.use[1])
    plot.data$plot.x <- object@meta.data[, idx]
  } else if (item.use[1] %in% colnames(object@dm@eigenvectors) ) {
    idx <- which(colnames(object@dm@eigenvectors) == item.use[1])
    plot.data$plot.x <- object@dm@eigenvectors[, idx]
  } else if (item.use[1] %in% colnames(object@pseudotime) ) {
    idx <- which(colnames(object@pseudotime) == item.use[1])
    plot.data$plot.x <- object@pseudotime[, idx]
  }

  # check paramter for plot.y
  if (item.use[2] %in% colnames(object@pca.load) ) {
    idx <- which(colnames(object@pca.load) == item.use[2])
    plot.data$plot.y <- object@pca.load[, idx]
  } else if (item.use[2] %in% colnames(object@tsne.value) ) {
    idx <- which(colnames(object@tsne.value) == item.use[2])
    plot.data$plot.y <- object@tsne.value[, idx]
  } else if (item.use[2] %in% colnames(object@gate.data) ) {
    idx <- which(colnames(object@gate.data) == item.use[2])
    plot.data$plot.y <- object@gate.data[, idx]
  } else if (item.use[2] %in% colnames(object@meta.data) ) {
    idx <- which(colnames(object@meta.data) == item.use[2])
    plot.data$plot.y <- object@meta.data[, idx]
  } else if (item.use[2] %in% colnames(object@dm@eigenvectors) ) {
    idx <- which(colnames(object@dm@eigenvectors) == item.use[2])
    plot.data$plot.y <- object@dm@eigenvectors[, idx]
  } else if (item.use[2] %in% colnames(object@pseudotime) ) {
    idx <- which(colnames(object@pseudotime) == item.use[2])
    plot.data$plot.y <- object@pseudotime[, idx]
  }

  # check paramter for color.by
  if (!color.by %in% colname.info ) {
    message(paste0(Sys.time(), " [WARNING] invalided color.by parameter "))
    plot.data$color.by <- 0
  } else {
    if (color.by %in% colnames(object@pca.load) ) {
      idx <- which(colnames(object@pca.load) == color.by)
      plot.data$color.by <- object@pca.load[, idx]
    } else if (color.by %in% colnames(object@tsne.value) ) {
      idx <- which(colnames(object@tsne.value) == color.by)
      plot.data$color.by <- object@tsne.value[, idx]
    } else if (color.by %in% colnames(object@gate.data) ) {
      idx <- which(colnames(object@gate.data) == color.by)
      plot.data$color.by <- object@gate.data[, idx]
    } else if (color.by %in% colnames(object@meta.data) ) {
      idx <- which(colnames(object@meta.data) == color.by)
      plot.data$color.by <- object@meta.data[, idx]
    } else if (color.by %in% colnames(object@dm@eigenvectors) ) {
      idx <- which(colnames(object@dm@eigenvectors) == color.by)
      plot.data$color.by <- object@dm@eigenvectors[, idx]
    } else if (color.by %in% colnames(object@pseudotime) ) {
      idx <- which(colnames(object@pseudotime) == color.by)
      plot.data$color.by <- object@pseudotime[, idx]
    } else {
      message(paste0(Sys.time(), " [WARNING] invalided color.by parameter "))
      plot.data$color.by <- 0
    }
  }

  if (is.numeric(plot.data$color.by)) {
    gg <- ggplot(plot.data, aes(x=plot.x, y=plot.y, colour = color.by)) + geom_point(size = size, alpha = alpha)
    gg <- gg + scale_colour_gradientn(colours = color.theme)
  } else {
    gg <- ggplot(plot.data, aes(x=plot.x, y=plot.y, colour = color.by)) + geom_point(size = size, alpha = alpha)
    gg <- gg + scale_color_manual(values = color.theme)
  }

  gg <- gg + theme_base()
  gg <- gg + labs(x = item.use[1], y = item.use[2], title = paste0(main))

  return(gg)
}



#'
#' Visualization of 2D data
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 3D plot, axes x and y and z must be numeric.
#' @param color.by character. Dot color by which character. It can be one of the column
#'     of raw.meta.data, or it can be just "density" (the default value).
#' @param size numeric. size of the dot
#' @param scale.y numeric. scale of y axis related to x- and z axis
#' @param angle numberic. angle between x and y axis (Attention: result depends on scaling).
#' @param main character. title of the plot
#' @param color.theme character. Library of color theme
#'
#' @import ggplot2 scatterplot3d
#'
#' @export
#'
#' @examples
#'
plot3D <- function(object,
                   item.use = c("PC1", "PC2", "PC3"),
                   color.by = "stage",
                   scale.y = 0.9,
                   angle = 40,
                   size = 1,
                   main = "3D plot of FSPY",
                   color.theme = c("#f2de00", "#8bd129", "#33CCFF", "#0066CC", "#CC66FF", "#FF3300")) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")

  # checking items
  colname.info <- c(colnames(object@pca.load),
                    colnames(object@tsne.value),
                    colnames(object@gate.data),
                    colnames(object@meta.data),
                    colnames(object@dm@eigenvectors),
                    colnames(object@pseudotime))

  if ( (sum(!item.use %in% colname.info) > 0) | (length(item.use) != 3) ) stop(paste0(Sys.time(), " [ERROR] invalided colname.info paramter "))

  plot.data <- data.frame(object@meta.data)

  # check paramter for plot.x
  if (item.use[1] %in% colnames(object@pca.load) ) {
    idx <- which(colnames(object@pca.load) == item.use[1])
    plot.data$plot.x <- object@pca.load[, idx]
  } else if (item.use[1] %in% colnames(object@tsne.value) ) {
    idx <- which(colnames(object@tsne.value) == item.use[1])
    plot.data$plot.x <- object@tsne.value[, idx]
  } else if (item.use[1] %in% colnames(object@gate.data) ) {
    idx <- which(colnames(object@gate.data) == item.use[1])
    plot.data$plot.x <- object@gate.data[, idx]
  } else if (item.use[1] %in% colnames(object@meta.data) ) {
    idx <- which(colnames(object@meta.data) == item.use[1])
    plot.data$plot.x <- object@meta.data[, idx]
  } else if (item.use[1] %in% colnames(object@dm@eigenvectors) ) {
    idx <- which(colnames(object@dm@eigenvectors) == item.use[1])
    plot.data$plot.x <- object@dm@eigenvectors[, idx]
  } else if (item.use[1] %in% colnames(object@pseudotime) ) {
    idx <- which(colnames(object@pseudotime) == item.use[1])
    plot.data$plot.x <- object@pseudotime[, idx]
  }

  # check paramter for plot.y
  if (item.use[2] %in% colnames(object@pca.load) ) {
    idx <- which(colnames(object@pca.load) == item.use[2])
    plot.data$plot.y <- object@pca.load[, idx]
  } else if (item.use[2] %in% colnames(object@tsne.value) ) {
    idx <- which(colnames(object@tsne.value) == item.use[2])
    plot.data$plot.y <- object@tsne.value[, idx]
  } else if (item.use[2] %in% colnames(object@gate.data) ) {
    idx <- which(colnames(object@gate.data) == item.use[2])
    plot.data$plot.y <- object@gate.data[, idx]
  } else if (item.use[2] %in% colnames(object@meta.data) ) {
    idx <- which(colnames(object@meta.data) == item.use[2])
    plot.data$plot.y <- object@meta.data[, idx]
  } else if (item.use[2] %in% colnames(object@dm@eigenvectors) ) {
    idx <- which(colnames(object@dm@eigenvectors) == item.use[2])
    plot.data$plot.y <- object@dm@eigenvectors[, idx]
  } else if (item.use[2] %in% colnames(object@pseudotime) ) {
    idx <- which(colnames(object@pseudotime) == item.use[2])
    plot.data$plot.y <- object@pseudotime[, idx]
  }

  # check paramter for plot.y
  if (item.use[3] %in% colnames(object@pca.load) ) {
    idx <- which(colnames(object@pca.load) == item.use[3])
    plot.data$plot.z <- object@pca.load[, idx]
  } else if (item.use[3] %in% colnames(object@tsne.value) ) {
    idx <- which(colnames(object@tsne.value) == item.use[3])
    plot.data$plot.z <- object@tsne.value[, idx]
  } else if (item.use[3] %in% colnames(object@gate.data) ) {
    idx <- which(colnames(object@gate.data) == item.use[3])
    plot.data$plot.z <- object@gate.data[, idx]
  } else if (item.use[3] %in% colnames(object@meta.data) ) {
    idx <- which(colnames(object@meta.data) == item.use[3])
    plot.data$plot.z <- object@meta.data[, idx]
  } else if (item.use[3] %in% colnames(object@dm@eigenvectors) ) {
    idx <- which(colnames(object@dm@eigenvectors) == item.use[3])
    plot.data$plot.z <- object@dm@eigenvectors[, idx]
  } else if (item.use[3] %in% colnames(object@pseudotime) ) {
    idx <- which(colnames(object@pseudotime) == item.use[3])
    plot.data$plot.z <- object@pseudotime[, idx]
  }

  # check paramter for color.by
  if (!color.by %in% colname.info ) {
    message(paste0(Sys.time(), " [WARNING] invalided color.by parameter "))
    plot.data$color.by <- 0
  } else {
    if (color.by %in% colnames(object@pca.load) ) {
      idx <- which(colnames(object@pca.load) == color.by)
      plot.data$color.by <- object@pca.load[, idx]
    } else if (color.by %in% colnames(object@tsne.value) ) {
      idx <- which(colnames(object@tsne.value) == color.by)
      plot.data$color.by <- object@tsne.value[, idx]
    } else if (color.by %in% colnames(object@gate.data) ) {
      idx <- which(colnames(object@gate.data) == color.by)
      plot.data$color.by <- object@gate.data[, idx]
    } else if (color.by %in% colnames(object@meta.data) ) {
      idx <- which(colnames(object@meta.data) == color.by)
      plot.data$color.by <- object@meta.data[, idx]
    } else if (color.by %in% colnames(object@dm@eigenvectors) ) {
      idx <- which(colnames(object@dm@eigenvectors) == color.by)
      plot.data$color.by <- object@dm@eigenvectors[, idx]
    } else if (color.by %in% colnames(object@pseudotime) ) {
      idx <- which(colnames(object@pseudotime) == color.by)
      plot.data$color.by <- object@pseudotime[, idx]
    } else {
      message(paste0(Sys.time(), " [WARNING] invalided color.by parameter "))
      plot.data$color.by <- 0
    }
  }

  if (is.numeric(plot.data$color.by)) {
    color.lib <- colorRampPalette(color.theme)(102)
    plot.data$color.by.sd <- plot.data$color.by - min(plot.data$color.by)
    plot.data$color.by.3d <- color.lib[ ceiling( plot.data$color.by.sd/max(plot.data$color.by.sd) * 100 ) + 1 ]
  } else {
    color.num <- unique(plot.data$color.by)
    color.lib <- colorRampPalette(color.theme)(length(color.num))
    plot.data$color.by.3d <- color.lib[match(plot.data$color.by, color.num)]
  }

  scatterplot3d(x = plot.data$plot.x, y = plot.data$plot.y, z = plot.data$plot.z,
                color = plot.data$color.by.3d,
                pch = 16, cex.symbols = size,
                scale.y = scale.y, angle = angle,
                xlab = item.use[1], ylab = item.use[2], zlab = item.use[3],
                main = main,
                col.axis = "#444444", col.grid = "#CCCCCC")

}

#'
#' plot SOM network
#'
#' @param object an FSPY object
#' @param cex.size numeric. size cex of the dot
#' @param color.by numeric. size color theme of the dot
#' @param size.by numeric. size theme of the dot
#' @param show.node.name logical. whether to show node name
#' @param color.theme character. Library of color theme
#'
#' @export
#'
#' @examples
#'
plotSOM <- function(object,
                    cex.size = 1,
                    color.by = "default",
                    size.by = "cell.percent",
                    show.node.name = F,
                    color.theme =  c("#0099FF", "#FFFF33", "#FF0033")) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  if (is.null(object@som)) stop(Sys.time(), " [ERROR] som is missing, please run runSOM first!")

  layout.attr <- object@som$layout
  node.attr <- object@som$node.attr
  edge.attr <- object@som$edge.attr


  if (!size.by %in% colnames(node.attr)) {
    message(Sys.time(), " [WARNING] invalid size.by attribute")
    size.by <- "cell.percent"
  }
  if (color.by != "default") {
    if (!color.by %in% colnames(node.attr)) {
      message(Sys.time(), " [WARNING] invalid size.by attribute")
      color.by <- "default"
      color <- 0
    } else {
      color <- node.attr[, which(colnames(node.attr) == color.by)]
    }
  } else {
    color <- 0
  }

  size <- node.attr[, which(colnames(node.attr) == size.by)]

  gg <- ggplot()
  gg <- gg + geom_segment(aes(x = edge.attr$from.x, y = edge.attr$from.y, xend = edge.attr$to.x, yend = edge.attr$to.y))
  gg <- gg + geom_point(mapping = aes(x = layout.attr$pos.x, y = layout.attr$pos.y, size = size, color = color))

  if (is.numeric(color)) {
    gg <- gg + scale_colour_gradientn(colours = color.theme)
  } else {
    gg <- gg + scale_color_manual(values = color.theme)
  }

  gg <- gg + scale_size(range = c(0, 6) * cex.size)

  if (show.node.name) gg <- gg + geom_text(aes(layout.attr$pos.x, y = layout.attr$pos.y, label = rownames(node.attr) ), check_overlap = TRUE, size = 3 * cex.size)
  gg <- gg + theme_void()
  gg <- gg + labs(x = "", y = "", title = paste0("SOM plot, color.by: ", color.by, ", size.by: ", size.by))

  return(gg)

}


#'
#' plot heatmap of FSPY
#'
#' @param object an FSPY object
#' @param color.by character.
#' @param main character. title of the plot
#' @param color.theme character. Library of color theme
#'
#'
#' @import pheatmap
#'
#' @export
#'
#' @examples
#'
plotPseudotimeDensity <- function(object, color.by = "stage",
                                  main = "Density of pseudotime",
                                  adjust = 0.5,
                                  color.theme = c("#f2de00", "#8bd129", "#33CCFF", "#0066CC", "#CC66FF", "#FF3300")) {
  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")

  # checking items
  colname.info <- colnames(object@meta.data)

  if (!color.by %in% colname.info) {
    stop(Sys.time(), " [ERROR] ", color.by, " is not in colnames of meta.data")
  } else {
    plot.data <- data.frame(
      cell = object@meta.data$cell,
      pseudotime = object@meta.data$pseudotime,
      color.by = object@meta.data[, which(colnames(object@meta.data) == color.by)]
    )
  }

  if (!is.factor(plot.data$color.by)) plot.data$color.by <- as.factor(as.character(plot.data$color.by))

  gg <- ggplot(plot.data, aes(x=pseudotime, colour = color.by))
  color.num <- unique(plot.data$color.by)
  color.lib <- colorRampPalette(color.theme)(length(color.num))
  gg <- gg + scale_color_manual(values = color.lib)
  gg <- gg + geom_density(adjust = adjust)

  gg <- gg + theme_base()
  gg <- gg + labs(title = paste0(main))

  return(gg)
}

#'
#' plot som tree of FSPY
#'
#' @param object an FSPY object
#' @param color.by character.
#' @param main character. title of the plot
#' @param color.theme character. Library of color theme
#'
#' @export
#'
#' @examples
#'
#'
#'
plotSOMtree <- function(object, root.id = NULL,
                        cex.size = 1,
                        color.by = "default",
                        size.by = "cell.percent",
                        show.node.name = F,
                        color.theme =  c("#f2de00", "#8bd129", "#33CCFF", "#0066CC", "#CC66FF", "#FF3300")) {
  if (is.null(root.id)) root.id <- which( object@som$node.attr$pseudotime == min( object@som$node.attr$pseudotime) )
  layout <- as.data.frame(igraph::layout_as_tree(object@som$graph, root = root.id ))

  layout.attr <- data.frame(pos.x = -1 * (layout$V2),
                            pos.y = layout$V1)
  node.attr <- object@som$node.attr
  edge.attr <- object@som$edge.attr

  edge.attr$from.x <- layout.attr$pos.x[match(edge.attr$from, rownames(layout))]
  edge.attr$from.y <- layout.attr$pos.y[match(edge.attr$from, rownames(layout))]
  edge.attr$to.x <- layout.attr$pos.x[match(edge.attr$to, rownames(layout))]
  edge.attr$to.y <- layout.attr$pos.y[match(edge.attr$to, rownames(layout))]


  if (!size.by %in% colnames(node.attr)) {
    message(Sys.time(), " [WARNING] invalid size.by attribute")
    size.by <- "cell.percent"
  }
  if (color.by != "default") {
    if (!color.by %in% colnames(node.attr)) {
      message(Sys.time(), " [WARNING] invalid size.by attribute")
      color.by <- "default"
      color <- 0
    } else {
      color <- node.attr[, which(colnames(node.attr) == color.by)]
    }
  } else {
    color <- 0
  }

  size <- node.attr[, which(colnames(node.attr) == size.by)]

  gg <- ggplot()
  gg <- gg + geom_segment(aes(x = edge.attr$from.x, y = edge.attr$from.y, xend = edge.attr$to.x, yend = edge.attr$to.y))
  gg <- gg + geom_point(mapping = aes(x = layout.attr$pos.x, y = layout.attr$pos.y, size = size, color = color))

  if (is.numeric(color)) {
    gg <- gg + scale_colour_gradientn(colours = color.theme)
  } else {
    gg <- gg + scale_color_manual(values = color.theme)
  }

  gg <- gg + scale_size(range = c(0, 6) * cex.size)

  if (show.node.name) gg <- gg + geom_text(aes(layout.attr$pos.x, y = layout.attr$pos.y, label = rownames(node.attr) ), check_overlap = TRUE, size = 3 * cex.size)
  gg <- gg + theme_void()
  gg <- gg + labs(x = "", y = "", title = paste0("SOM plot, color.by: ", color.by, ", size.by: ", size.by))

  return(gg)
}











