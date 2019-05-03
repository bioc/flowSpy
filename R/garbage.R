
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






  #' Calculate pseudotime by 'flooding'
  #'
  #' This calculates pseudotime by performing a probabilistic breadth-first search
  #' of the k-nearest neighbor graph that was used to generate the diffusion map
  #' calculated on the data. The results of this function should be passed to
  #' \code{\link{floodPseudotimeProcess}} to convert them into pseudotime.
  #'
  #' @keywords internal
  #'
  #' @references
  #'   This code is strongly based on the \code{\link[URD]{floodPseudotime}} function,
  #'   which is developed by Jeffrey A. Farrell. See paper Jeffrey A. Farrell, et. al.,
  #'   Single-cell reconstruction of developmental trajectories during zebrafish
  #'   embryogenesis, Science, 2018.
  #'
  #'
  pseudotimeFlood <- function(object, n = 10, minimum.cells.flooded = 2, thread = 1, verbose = F) {

    if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
    if (is.null(object@dm)) stop(Sys.time(), " [ERROR] dm in FSPY is missing, please run runDiffusionMap first")

    tm.flood <- object@dm@transitions
    tm.sum.max <- max(apply(tm.flood, 1, sum))
    tm.flood <- tm.flood / tm.sum.max


    floods <- as.data.frame(lapply(1:n, function(x) {
      if (verbose) print(paste("Starting flood number", x))
      floodPseudotimeCalc(object, minimum.cells.flooded = minimum.cells.flooded, tm.flood = tm.flood)
    }))

    colnames(floods) <- 1:n

    return(floods)
  }


  #' A single iteration of flood pseudotime calculation.
  #'
  #' @param object An FSPY object
  #' @param minimum.cells.flooded (Numeric) Stop pseudotime calculation when fewer cells are newly visited in a given iteration, and assign remaining unvisited cells pseudotime NA.
  #' @param tm.flood (Sparse or full matrix) A sparse matrix of normalized transition probabilities. If unprovided (\code{NULL}), this will be calculated automatically.
  #'
  #' @return (Numeric vector) The iteration that newly visited each cell.
  #'
  #' @keywords internal
  #'
  #' @references
  #'   This code is strongly based on the \code{\link[URD]{floodPseudotime}} function,
  #'   which is developed by Jeffrey A. Farrell. See paper Jeffrey A. Farrell, et. al.,
  #'   Single-cell reconstruction of developmental trajectories during zebrafish
  #'   embryogenesis, Science, 2018.
  #'
  #'
  floodPseudotimeCalc <- function(object, minimum.cells.flooded, tm.flood) {
    pseudotime <- rep(NA, dim(tm.flood)[1])
    names(pseudotime) <- rownames(tm.flood)
    i <- 0

    # Initialize with starting cells
    cells.visited <- object@root.cells
    cells.not.visited <- setdiff(names(pseudotime), cells.visited)
    pseudotime[root.cell] <- 0

    total.cells <- dim(tm.flood)[1]
    cells.newly.visited <- rep(NA, minimum.cells.flooded)

    # Loop until all cells have a pseudotime
    while (length(cells.visited) < (total.cells-1) & length(cells.newly.visited) >= minimum.cells.flooded) {
      # Increment step counter
      i <- i + 1
      # Calculate visitation probability for each cell
      visitation.prob <- apply(tm.flood[cells.not.visited, cells.visited], 1, combine.probs)
      # Figure out which cells are newly visited
      cells.newly.visited <- names(visitation.prob)[which(rbinom(length(visitation.prob), 1, visitation.prob) == 1)]
      # Record the visited cells
      cells.visited <- c(cells.visited, cells.newly.visited)
      cells.not.visited <- setdiff(cells.not.visited, cells.newly.visited)
      pseudotime[cells.newly.visited] <- i
    }

    return(pseudotime)
  }

  #' Combine probabilities
  #'
  #' Calculates cumulative probability
  #'
  #' @param x (Numeric Vector) Probabilities
  #'
  #' @return (Numeric) Cumulative probability
  #'
  #' @keywords internal
  #'
  #'
  combine.probs <- function(x) {
    1-prod(1-x)
  }

  #' Process Flood Pseudotime
  #'
  #' This processes the values returned by \code{\link{floodPseudotime}} and stores
  #' the result as a pseudotime in an URD object.
  #'
  #' @export
  #'
  #' @examples
  #'
  #' @references
  #'   This code is strongly based on the \code{\link[URD]{floodPseudotime}} function,
  #'   which is developed by Jeffrey A. Farrell. See paper Jeffrey A. Farrell, et. al.,
  #'   Single-cell reconstruction of developmental trajectories during zebrafish
  #'   embryogenesis, Science, 2018.
  #'
  pseudotimeProcess <- function(object, pseudotime.fun = mean, stability.div = 10,
                                n = 10, minimum.cells.flooded = 2, thread = 1, verbose = F) {

    floods <- pseudotimeFlood(object, n = n, minimum.cells.flooded = minimum.cells.flooded, thread = thread, verbose = verbose)

    # Turn flood positions into relative positions
    # (e.g. instead of ordinal positioning, fractional position from 0 to 1)
    floods <- sweep(floods, 2, apply(floods, 2, max, na.rm=T), "/")
    # Divide the floods into sections for stability calculations
    floods.in.division <- ceiling(1:stability.div / stability.div * dim(floods)[2])
    ## CALCULATE PSEUDOTIME AND VISIT FREQUENCY WITH INCREASING AMOUNTS OF DATA TO TEST WHEN IT STABILIZES
    walks.per.cell <- as.data.frame(lapply(floods.in.division, function(n.floods) {
      if (n.floods == 1) return(as.numeric(!is.na(floods[,1])))
      else return (apply(floods[,1:n.floods], 1, function(x) sum(as.numeric(!is.na(x)))))
    }))
    pseudotime.stability <- as.data.frame(lapply(floods.in.division, function(n.floods) {
      if (n.floods == 1) return(floods[,1])
      else return (apply(floods[,1:n.floods], 1, pseudotime.fun, na.rm=T))
    }))
    colnames(walks.per.cell) <- seq(1, ncol(floods), length.out=stability.div)
    colnames(pseudotime.stability) <- seq(1, ncol(floods), length.out=stability.div)

    for (i in 1:n) pseudotime.stability[which(is.na(pseudotime.stability[, i])), i] = 0
    pseudotime <- data.frame(
      pseudotime.turn = pseudotime.stability
    )

    # Store the results: after calculating with all walks, visit frequencies in diff.data
    final.visit.freq <- walks.per.cell[rownames(floods),dim(walks.per.cell)[2]]

    pseudotime$visitfreq.raw.pseudotime <- final.visit.freq
    pseudotime$visitfreq.log.pseudotime <- log10(final.visit.freq + 1)

    pseudotime$pseudotime <- rowMeans(pseudotime.stability, na.rm = T)
    min.cutoff <- min(pseudotime$pseudotime[pseudotime$pseudotime > 0], na.rm = T)
    pseudotime$pseudotime[pseudotime$pseudotime > 0] <- pseudotime$pseudotime[pseudotime$pseudotime > 0] - min.cutoff
    pseudotime$pseudotime <- (pseudotime$pseudotime / max(pseudotime$pseudotime)) ** 0.5

    object@meta.data$pseudotime <- pseudotime$pseudotime
    object@som$node.attr$pseudotime <- unlist(lapply(rownames(object@som$node.attr), function(x) median( pseudotime$pseudotime[which(object@meta.data$som.node.id == x)] )))

    return(object)
  }





}
