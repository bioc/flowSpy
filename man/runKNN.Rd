% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn.R
\name{runKNN}
\alias{runKNN}
\title{Calculate k-nearest neighbors of FSPY}
\usage{
runKNN(
  object,
  given.mat = NULL,
  knn = 30,
  knn.replace = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{object}{an FSPY object}

\item{given.mat}{matrix. Given matrix to run knn}

\item{knn}{numeric. Number of k-nearest neighbors.}

\item{knn.replace}{logic. Whether to replace knn in FSPY object}

\item{verbose}{logical. Whether to print calculation progress.}

\item{...}{Parameters passing to \code{\link[BiocNeighbors]{findKNN}} function}
}
\value{
An FSPY object with knn, knn.index and knn.distance information.
}
\description{
Calculates and stores a k-nearest neighbor graph based on Euclidean
   distance with (KMKNN) algorithm using log-transformed signaling matrix of
   flow cytometry data. The base function are base on \code{\link[BiocNeighbors]{findKNN}}.
}
\examples{
if (FALSE) {

fspy <- runKNN(fspy)

}


}
\seealso{
\code{\link[BiocNeighbors]{findKNN}}
}
