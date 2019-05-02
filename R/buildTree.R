#
#
# build SOM network
#
buildTree <- function(object, method = "euclidean") {

  som.net <- list()
  adjacency <- stats::dist(flowsom$codes, method = method)
  fullGraph <- igraph::graph.adjacency(as.matrix(adjacency),
                                       mode = "undirected",
                                       weighted = TRUE)
  fullGraph <- simplify(fullGraph)
  som.net$graph <- igraph::minimum.spanning.tree(fullGraph)
  som.net$layout <- as.data.frame(igraph::layout.kamada.kawai(som.net$graph))
  colnames(som.net$layout) <- c("Pos.x", "Pos.y")

  som.net$node.attr <- data.frame(cell.num = as.vector(table(flowsom$mapping[, 1])),
                                  cell.percent = as.vector(table(flowsom$mapping[, 1])/dim(flowsom$mapping)[1]),
                                  aggregate(object@log.data, list(som.id = object@meta.data$som.id), mean))

  idx <- match(c("som.id", "stage"), colnames(object@meta.data))

  cell.percent <- matrix(table(object@meta.data[, idx]), nrow = max(object@meta.data$som.id))
  colnames(cell.percent) <- levels(object@meta.data$stage)
  cell.percent.stage <- cell.percent / som.net$node.attr$cell.num
  colnames(cell.percent.stage) <- paste0(levels(object@meta.data$stage), ".som.percent")

  som.net$node.attr <- cbind(som.net$node.attr, cell.percent, cell.percent.stage)

  som.net$edge.attr <- igraph::as_data_frame(som.net$graph, what="edges")
  som.net$edge.attr$from.x <- som.net$layout$Pos.x[match(som.net$edge.attr$from, rownames(som.net$layout))]
  som.net$edge.attr$from.y <- som.net$layout$Pos.y[match(som.net$edge.attr$from, rownames(som.net$layout))]
  som.net$edge.attr$to.x <- som.net$layout$Pos.x[match(som.net$edge.attr$to, rownames(som.net$layout))]
  som.net$edge.attr$to.y <- som.net$layout$Pos.y[match(som.net$edge.attr$to, rownames(som.net$layout))]

  return(som.net)
}











