#'
#' buildTree
#'
#'
buildTree <- function(object, method = "euclidean",
                      cluster.type = "som",
                      dim.type = "tsne", dim.use = 1:2,
                      verbose = T) {

  if (verbose) message(Sys.time(), " [INFO] Calculating buildTree.")
  if (missing(object)) stop(Sys.time(), " [ERROR] FSPY object is missing.")
  if (dim.type %in% c("tsne", "tSNE", "TSNE", "t-SNE","t_SNE", "t") ) {
    dim.name <- paste0("tSNE", dim.use)
    tree.mat <- object@tsne.value[, dim.name]
  } else if ( dim.type %in% c("PCA", "pca", "p") ) {
    dim.name <- paste0("PC", dim.use)
    tree.mat <- object@pca.value[, dim.name]
  } else if (dim.type %in% c("dc", "diffusionmap", "diffusion-map", "destiny", "d")) {
    dim.name <- paste0("DC", dim.use)
    tree.mat <- object@dm@eigenvectors[, dim.name]
  } else if (dim.type %in% c("umap", "UMAP", "u")) {
    dim.name <- paste0("UMAP", dim.use)
    tree.mat <- object@umap.value[, dim.name]
  } else {
    if (verbose) message(Sys.time(), " [INFO] the log data will be used to calculate trajectory")
    tree.mat <- object@log.data
  }

  cluster.name <- paste0(cluster.type, ".id")
  if (! cluster.name %in% colnames(object@meta.data)) {
    stop(Sys.time(), " [ERROR] Invalid cluster.name ")
  }
  mst.mat <- aggregate(tree.mat, list(cluster = object@meta.data[, cluster.name]), mean)
  rownames(mst.mat) <- mst.mat[, 1]

  adjacency <- stats::dist(mst.mat[, -1], method = method)
  fullGraph <- igraph::graph.adjacency(as.matrix(adjacency),
                                       mode = "undirected",
                                       weighted = TRUE)
  fullGraph <- simplify(fullGraph)
  tree.graph <- igraph::minimum.spanning.tree(fullGraph)

  object@network <- list(mst = tree.graph,
                         dist = adjacency,
                         method = method,
                         cluster.type = cluster.type,
                         dim.type = dim.type,
                         dim.use = dim.use,
                         tree.mat = tree.mat,
                         mst.mat = mst.mat)

  if (verbose) message(Sys.time(), " [INFO] Calculating buildTree completed.")
  return(object)
}











