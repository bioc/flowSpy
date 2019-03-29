#'
#' definition of root cells
#'
#' @param object an FSPY object
#' @param root.cells vector, cell name of the root cells
#'
#' @export
#'
#' @examples
#'
#'
#'
defRootCells <- function(object, root.cells = NULL) {
  if (length(object@root.cells) != 0) warning(Sys.time(), " [WARNING] root cells in FSPY object exist, they will be replaced.")

  if (!is.vector(root.cells)) stop(Sys.time(), " [ERROR] root.cell must be a vector")

  root.cells <- root.cells[root.cells %in% object@meta.data$cell]
  if ( length(root.cells) == 0 ) {
    stop(Sys.time(), " [ERROR] root.cell not in meta.data")
  } else {
    object@root.cells <- root.cells
  }

  return(object)
}


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














