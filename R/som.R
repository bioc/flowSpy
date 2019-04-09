#'
#' calculation SOM in FSPY object
#'
#' @description
#' Build a self-organizing map
#'
#' @param object  an FSPY object
#' @param xdim  Width of the grid.
#' @param ydim  Hight of the grid.
#' @param rlen  Number of times to loop over the training data for each MST
#' @param mst   Number of times to build an MST
#' @param alpha Start and end learning rate
#' @param radius Start and end radius
#' @param init  Initialize cluster centers in a non-random way
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev,
#'              4=cosine)
#' @param codes Cluster centers to start with
#' @param importance array with numeric values. Parameters will be scaled
#'                   according to importance
#' @param method the distance measure to be used. This must be one of "euclidean",
#'      "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'      Any unambiguous substring can be given. See \code{\link[stats]{dist}}
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return an FSPY object
#' @seealso \code{\link{BuildSOM}}
#'
#' @references This code is strongly based on the \code{\link[FlowSOM]{SOM}} function.
#'             Which is developed by Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2018).
#'
#' @import FlowSOM
#' @importFrom igraph graph.adjacency minimum.spanning.tree layout.kamada.kawai
#'
#' @export
#'
#' @examples
#'
runSOM <- function(object, xdim = 5, ydim = 5, rlen = 8, mst = 1,
                   alpha = c(0.05,  0.01), radius = 1, init = FALSE,
                   distf = 2, codes = NULL, importance = NULL,
                   method = "euclidean", verbose= T) {

  if (verbose) message(Sys.time(), " [INFO] Calculating FlowSOM.")
  # flowSOM
  flowset <- as.matrix(object@log.data)
  flowsom <- FlowSOM::SOM(flowset,
                          xdim = xdim, ydim = ydim, rlen = rlen, mst = mst,
                          alpha = alpha[1], radius = radius,
                          init = init,
                          distf = distf, silent = verbose,
                          codes = codes, importance = importance)


  # generating som network
  object@meta.data$som.id <- flowsom$mapping[, 1]
  object@meta.data$som.value <- flowsom$mapping[, 2]
  object@som <- flowsom

  object@som.network <- buildSOMnet(flowsom, object, method = method)

  if (verbose) message(Sys.time(), " [INFO] Calculating FlowSOM completed.")
  return(object)
}


#
# build SOM network
#
buildSOMnet <- function(flowsom, object, method = "euclidean") {

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











