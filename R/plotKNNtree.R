#'
#' plot Minimum spanning tree from knn cluster of FSPY
#'
#' @param object an FSPY object
#' @param cex.size numeric. size cex of the dot
#' @param color.by numeric. size color theme of the dot
#' @param size.by numeric. size theme of the dot
#' @param show.node.name logical. whether to show node name
#' @param color.theme character. Library of color theme
#'
#'
plotKNNtree <- function(object,
                        cex.size = 1,
                        color.by = "default",
                        size.by = "cell.percent",
                        layout = "mst",
                        show.node.name = F,
                        color.theme =  NULL) {
  if (missing(object)) {
    stop(Sys.time(), "[ERROR] object is missing")
  }
  if (is.null(object@trunk.network$trunk.graph)) {
    stop(Sys.time(), "[ERROR] trunk is missing in object")
    layout <- as.data.frame(igraph::layout_with_fr(object@trunk.network$trunk.spanning.tree))
    colnames(layout) <- paste0("Pos", 1:dim(layout)[2])
  }

  if (!size.by %in% colnames(object@trunk.network$trunk.marker)) {
    warning(Sys.time(), " [WARNING] invalid size.by attribute")
    size.by <- "CD34"
  }
  if (color.by != "default") {
    if (!color.by %in% colnames(object@trunk.network$trunk.marker)) {
      warning(Sys.time(), " [WARNING] invalid size.by attribute")
      color.by <- "default"
      color <- 0
    } else {
      color <- node.attr[, which(colnames(object@trunk.network$trunk.marker) == color.by)]
    }
  } else {
    color <- 0
  }


  edge.attr <- as.data.frame(as_edgelist(object@trunk.network$trunk.spanning.tree))
  colnames(edge.attr) <- c("from", "to")

  edge.attr$from.x <- layout$Pos1[match(edge.attr$from, rownames(layout))]
  edge.attr$from.y <- layout$Pos2[match(edge.attr$from, rownames(layout))]
  edge.attr$to.x <- layout$Pos1[match(edge.attr$to, rownames(layout))]
  edge.attr$to.y <- layout$Pos2[match(edge.attr$to, rownames(layout))]

  gg <- ggplot()
  gg <- gg + geom_segment(aes(x = edge.attr$from.x, y = edge.attr$from.y, xend = edge.attr$to.x, yend = edge.attr$to.y))
  gg <- gg + geom_point(mapping = aes(x = layout$Pos1, y = layout$Pos2, size = size, color = color))

  if (is.numeric(color)) {
    gg <- gg + scale_colour_gradientn(colours = rainbow(10))
  } else {
    gg <- gg + scale_color_manual(values = rainbow(length(unique(color))))
  }

  gg <- gg + scale_size(range = c(0, 6) * cex.size)

  if (show.node.name) gg <- gg + geom_text(aes(layout$Pos1, y = layout$Pos2, label = rownames(layout) ), check_overlap = TRUE, size = 3 * cex.size)
  gg <- gg + theme_void()
  gg <- gg + labs(x = "", y = "", title = paste0("KNN tree plot, color.by: ", color.by, ", size.by: ", size.by))

  return(gg)

}





