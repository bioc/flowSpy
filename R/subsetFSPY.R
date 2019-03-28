#'
#' subset FSPY object
#'
#' @name subsetFSPY
#'
#' @description This subsets an FSPY object by given a list of cells or cluster id.
#'     This function will subset all results without recalculating them, such as knn,
#'     PCA, tSNE, umap and pseudotime. However, this function will worked after running
#'     \code{\link[flowSpy]{runKNN}}. For instance, you can choose recalculate PCA and
#'     tSNE and destiny scores by paramter recalculate.
#'
#' @param object An FSPY object
#' @param cells vector, Names of the cells to retain.
#' @param id.keep numeric or character. Name of the id to retain
#' @param id.name character. Name of the column which id.keep is in it.
#' @param recalculate logic. Whether to recalculate the PCA, tSNE, destiny and umap score.
#' @param knn numeric. When recalculate is TRUE.
#' @param verbose logic. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @examples
#'
#' @export
#'
subsetFSPY <- function(object, cells = NULL,
                       id.keep = NULL, id.name = NULL,
                       recalculate = F, knn = NA,
                       verbose = T) {
  # Make sure all cells are actually in the object
  cells.keep <- intersect(cells, rownames(object@log.data))

  object@log.data <- object@log.data[cells.keep, ]
  object@meta.data <- object@meta.data[cells.keep, ]

  if (!recalculate) {
    if(!any(dim(object@pca.scores) == 0))  object@pca.scores <- object@pca.scores[cells.keep, ]
    if(!any(dim(object@tsne.value) == 0))  object@tsne.value <- object@tsne.value[cells.keep, ]
  }

  return(object)
}



