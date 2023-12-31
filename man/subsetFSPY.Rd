% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetFSPY.R
\name{subsetFSPY}
\alias{subsetFSPY}
\title{subset FSPY object}
\usage{
subsetFSPY(object, cells = NULL, knn = NA, verbose = FALSE)
}
\arguments{
\item{object}{An FSPY object}

\item{cells}{vector, Names of the cells to retain.}

\item{knn}{numeric. If is NA, the KNN will be equal to the knn number in the input FSPY object.}

\item{verbose}{logic. Whether to print calculation progress.}
}
\value{
An FSPY object
}
\description{
This subsets an FSPY object by given a list of cells or cluster id.
    This function will subset all results without recalculating them, such as knn,
    PCA, tSNE, umap and pseudotime. For instance, you can choose recalculate PCA and
    tSNE and destiny scores by paramter recalculate.
}
\examples{

if (FALSE) {

cells <- test.meta.data$cell[which(test.meta.data$stage == "D0")]
sub.fspy <- subsetFSPY(fspy, cells = cells)
sub.fspy

}

}
