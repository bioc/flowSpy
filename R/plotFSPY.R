#'
#' plot heatmap of FSPY
#'
#' @param object an FSPY object
#' @param color.by character.
#' @param main character. title of the plot
#' @param color.theme character. Library of color theme
#'
#'
#' @import pheatmap
#'
#' @export
#'
#' @examples
#'
plotPseudotimeDensity <- function(object, color.by = "stage",
                                  main = "Density of pseudotime",
                                  adjust = 0.5,
                                  color.theme = NULL) {
  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")

  # checking items
  if ( !all("pseudotime" %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] pseudotime is not in plot.meta of FSPY, please run Pseudotime first.")

  if ( !all(color.by %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  item.use.idx <- match("pseudotime", colnames(object@plot.meta))
  color.by.idx <- match(color.by, colnames(object@plot.meta))

  plot.data <- data.frame(plot.x = object@plot.meta[, item.use.idx[1]],
                          color.by = object@plot.meta[, color.by.idx])


  if (!is.factor(plot.data$color.by)) plot.data$color.by <- as.factor(as.character(plot.data$color.by))

  gg <- ggplot(plot.data, aes(x=pseudotime, colour = color.by))
  color.lib <- rainbow(length(unique(plot.data$color.by)))
  gg <- gg + scale_color_manual(values = color.lib)
  gg <- gg + geom_density(adjust = adjust)

  gg <- gg + theme_base()
  gg <- gg + labs(title = paste0(main))

  return(gg)
}









