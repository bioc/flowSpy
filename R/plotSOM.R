#'
#' plot SOM network
#'
#' @name plotSOM
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
                    color.theme = NULL) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  if (is.null(object@som)) stop(Sys.time(), " [ERROR] som is missing, please run runSOM first!")

  layout.attr <- object@som$layout
  node.attr <- object@som$node.attr
  edge.attr <- object@som$edge.attr


  if (!size.by %in% colnames(node.attr)) {
    warning(Sys.time(), " [WARNING] invalid size.by attribute")
    size.by <- "cell.percent"
  }
  if (color.by != "default") {
    if (!color.by %in% colnames(node.attr)) {
      warning(Sys.time(), " [WARNING] invalid size.by attribute")
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
    gg <- gg + scale_colour_gradientn(colours = rainbow(10))
  } else {
    gg <- gg + scale_color_manual(values = rainbow(length(unique(color))))
  }

  gg <- gg + scale_size(range = c(0, 6) * cex.size)

  if (show.node.name) gg <- gg + geom_text(aes(x = layout.attr$pos.x, y = layout.attr$pos.y, label = rownames(node.attr) ), check_overlap = TRUE, size = 3 * cex.size)
  gg <- gg + theme_void()
  gg <- gg + labs(x = "", y = "", title = paste0("SOM plot, color.by: ", color.by, ", size.by: ", size.by))

  return(gg)

}




#'
#' plot som tree of FSPY
#'
#' @name plotSOMtree
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
                        color.theme =  NULL) {
  if (is.null(root.id)) {
    if ("pseudotime" %in% colnames(object@som$node.attr)) {
      root.id <- which( object@som$node.attr$pseudotime == min( object@som$node.attr$pseudotime) )
    } else {
      root.id <- 1
    }
  }
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
    warning(Sys.time(), " [WARNING] invalid size.by attribute")
    size.by <- "cell.percent"
  }
  if (color.by != "default") {
    if (!color.by %in% colnames(node.attr)) {
      warning(Sys.time(), " [WARNING] invalid size.by attribute")
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
    gg <- gg + scale_colour_gradientn(colours = rainbow(10))
  } else {
    gg <- gg + scale_color_manual(values = rainbow(length(unique(color))))
  }

  gg <- gg + scale_size(range = c(0, 6) * cex.size)

  if (show.node.name) gg <- gg + geom_text(aes(layout.attr$pos.x, y = layout.attr$pos.y, label = rownames(node.attr) ), check_overlap = TRUE, size = 3 * cex.size)
  gg <- gg + theme_void()
  gg <- gg + labs(x = "", y = "", title = paste0("SOM plot, color.by: ", color.by, ", size.by: ", size.by))

  return(gg)
}








