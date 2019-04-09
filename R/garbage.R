
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



}
