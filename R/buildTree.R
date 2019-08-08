#'
#' buildTree
#'
#' @name buildTree
#'
#' @param object an FSPY object
#' @param method character. Mehtod to build MST.
#' @param dim.type character.
#' @param dim.use numeric
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
#' @return An FSPY object
#'
#' @examples
#'
#' # build minimum spanning tree (MST) based on tsne
#' fspy <- buildTree(fspy, dim.type = "tsne", dim.use = 1:2, verbose = T)
#'
#' # Using PCA
#' fspy <- buildTree(fspy, dim.type = "pca", dim.use = 1:4, verbose = T)
#'
#' # Using UMAP
#' fspy <- buildTree(fspy, dim.type = "umap", dim.use = 1:2, verbose = T)
#'
#' # Using Diffusion Maps
#' fspy <- buildTree(fspy, dim.type = "dc", dim.use = 1:3, verbose = T)
#'
#'
buildTree <- function(object, method = "euclidean",
                      dim.type = "umap", dim.use = 1:2,
                      verbose = F) {

  if (verbose) message(Sys.time(), " [INFO] Calculating buildTree.")
  if (missing(object)) stop(Sys.time(), " [ERROR] FSPY object is missing.")
  if (dim.type %in% c("tsne", "tSNE", "TSNE", "t-SNE","t_SNE", "t") ) {
    dim.name <- paste0("tSNE_", dim.use)
    tree.mat <- object@tsne.value[, dim.name]
  } else if ( dim.type %in% c("PCA", "pca", "p") ) {
    dim.name <- paste0("PC_", dim.use)
    tree.mat <- object@pca.value[, dim.name]
  } else if (dim.type %in% c("dc", "diffusionmap", "diffusion-map", "destiny", "d")) {
    dim.name <- paste0("DC_", dim.use)
    tree.mat <- object@dm@eigenvectors[, dim.name]
  } else if (dim.type %in% c("umap", "UMAP", "u")) {
    dim.name <- paste0("UMAP_", dim.use)
    tree.mat <- object@umap.value[, dim.name]
  } else {
    if (verbose) message(Sys.time(), " [INFO] The log data will be used to calculate trajectory")
    tree.mat <- object@log.data[which(object@meta.data$dowsample == 1), ]
  }

  if (! "cluster.id" %in% colnames(object@meta.data)) {
    stop(Sys.time(), " [ERROR] Invalid cluster.id, please run runCluster first")
  }
  cluster.info <- object@meta.data$cluster.id[which(object@meta.data$dowsample == 1)]
  mst.mat <- aggregate(tree.mat, list(cluster = cluster.info), mean)
  rownames(mst.mat) <- mst.mat[, 1]

  adjacency <- stats::dist(mst.mat[, -1], method = method)
  fullGraph <- igraph::graph.adjacency(as.matrix(adjacency),
                                       mode = "undirected",
                                       weighted = TRUE)
  fullGraph <- simplify(fullGraph)
  tree.graph <- igraph::minimum.spanning.tree(fullGraph)

  # Storing network information
  object@network <- list(mst = tree.graph,
                         method = method,
                         dim.type = dim.type,
                         dim.use = dim.use,
                         mst.mat = mst.mat)

  # Initialization for root.cells and leaf cells
  if (verbose) message(Sys.time(), " [INFO] Initialization for root.cells and leaf cells")
  object@meta.data$is.root.cells <- 0
  object@meta.data$is.leaf.cells <- 0

  # update tree meta information
  object <- updateClustMeta(object, verbose = F)

  if (verbose) message(Sys.time(), " [INFO] Calculating buildTree completed.")
  return(object)
}











