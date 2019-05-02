
garbage <- function(object) {

  ############ KNN cluster
  if (F) {
    if (verbose) message(Sys.time(), " [INFO] start running knn.cluster. It will take some minutes if the dataset is too large. " )
    mat <- t(object@log.data)
    adj <- matrix(0, ncol(mat), ncol(mat))
    rownames(adj) <- colnames(adj) <- colnames(mat)
    for(i in seq_len(ncol(mat))) {
      adj[i, colnames(mat)[object@knn.index[i,]]] <- 1
    }
    g <- igraph::graph.adjacency(adj, mode="undirected")
    # remove self loops
    g <- simplify(g)
    ## identify communities
    km <- igraph::cluster_walktrap(g)
    # generation of trunk network
    object@network <- list(knn.G = g, knn.walktrap = km, adj = adj)
    object@meta.data$trunk.id <- km$membership

    if (verbose) message(Sys.time(), " [INFO] Add trunk ")
    object <- addTrunk(object)
    #if (verbose) message(Sys.time(), " [INFO] Add branch ")
    #object <- addBranch(object)

  } else {
    if (!"cluster.id" %in% colnames(object@meta.data)) {
      object@meta.data$trunk.id <- 0
    }
  }



  #
  # root.cells
  #
  addTrunk <- function(object) {
    trunk.network <- list()
    trunk.network$trunk.marker <- aggregate(object@log.data, list(cluster = object@meta.data$trunk.id), mean)
    rownames(trunk.network$trunk.marker) <- trunk.network$trunk.marker$cluster
    trunk.network$trunk.marker <- trunk.network$trunk.marker[, -1]
    trunk.network$trunk.dist <- stats::dist(trunk.network$trunk.marker, method = "euclidean")
    trunk.network$trunk.graph <- igraph::graph.adjacency(as.matrix(trunk.network$trunk.dist),
                                                         mode = "undirected",
                                                         weighted = TRUE)
    trunk.network$trunk.spanning.tree <- igraph::minimum.spanning.tree(trunk.network$trunk.graph)
    object@trunk.network <- trunk.network
    return(object)
  }


  #
  # root.cells
  #
  addBranch <- function(object) {
    trunk.id <- table(object@meta.data$trunk.id)
    object@meta.data$branch.id <- 0
    branch.sub.graph <- vector("list", length(names(trunk.id)))
    names(branch.sub.graph) <- names(trunk.id)
    for (trunk.id.sub in names(trunk.id)) {
      cells.sub.idx <- which(object@meta.data$trunk.id == trunk.id.sub)
      sub.i <- igraph::graph.adjacency(object@network$knn.G[cells.sub.idx, cells.sub.idx], mode="undirected")
      sub.i <- simplify(sub.i)
      branch.sub.graph[[trunk.id.sub]] <- sub.i
      km.sub <- igraph::cluster_walktrap(sub.i)
      object@meta.data$branch.id[cells.sub.idx] <- paste0(trunk.id.sub, "-" , km.sub$membership)
    }
    branch.network <- list(branch.sub.graph = branch.sub.graph)
    branch.network$branch.marker <- aggregate(object@log.data, list(cluster = object@meta.data$branch.id), mean)
    rownames(branch.network$branch.marker) <- branch.network$branch.marker$cluster
    branch.network$branch.marker <- branch.network$branch.marker[, -1]
    branch.network$branch.dist <- stats::dist(branch.network$branch.marker, method = "euclidean")
    branch.network$branch.graph <- igraph::graph.adjacency(as.matrix(branch.network$branch.dist),
                                                           mode = "undirected",
                                                           weighted = TRUE)
    branch.network$branch.spanning.tree <- igraph::minimum.spanning.tree(branch.network$branch.graph)
    object@branch.network <- branch.network
    return(object)
  }

  #'
  #' mclust
  #'
  #' @import mclust
  #'
  runMclust <- function(object) {
    mc.id <- Mclust(object@log.data)$classification

    object@meta.data$mc.id <- mc.id

    return(object)
  }



  tree.layout <- as.data.frame(igraph::layout.kamada.kawai(tree.graph))

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




}
