% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudotime.R
\name{runPseudotime}
\alias{runPseudotime}
\title{Calculation of Pseudotime}
\usage{
runPseudotime(
  object,
  mode = "undirected",
  dim.type = c("raw", "pca", "tsne", "dc", "umap"),
  dim.use = 1:2,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{object}{An FSPY object}

\item{mode}{character. Specifies how igraph should interpret the supplied matrix.
Possible values are: directed, undirected, upper, lower, max, min, plus.}

\item{dim.type}{character. Type of dimensionality reduction method used to calculate
pseudotime: raw, umap, tsne, dc and pca. By default is raw.}

\item{dim.use}{numeric. Dimensions used to calculate pseudotime}

\item{verbose}{logical. Whether to print calculation progress.}

\item{...}{Parameters passing to calculation function.}
}
\value{
An FSPY object
}
\description{
calculation of Pseudotime based on KNN
}
\examples{

if (FALSE) {

fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "raw")
fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "umap", dim.use = 1:2)
fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "tsne", dim.use = 1:2)
fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "dc", dim.use = 1:3)
fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "pca", dim.use = 1:3)

# tSNE plot colored by pseudotime
plot2D(fspy, item.use = c("tSNE_1", "tSNE_2"), category = "numeric",
       size = 1, color.by = "pseudotime") +
       scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))
# UMAP plot colored by pseudotime
plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), category = "numeric",
       size = 1, color.by = "pseudotime") +
       scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))
}

}
