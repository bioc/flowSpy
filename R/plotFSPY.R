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
                                  color.theme = c("#f2de00", "#8bd129", "#33CCFF", "#0066CC", "#CC66FF", "#FF3300")) {
  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")

  # checking items
  colname.info <- colnames(object@meta.data)

  if (!color.by %in% colname.info) {
    stop(Sys.time(), " [ERROR] ", color.by, " is not in colnames of meta.data")
  } else {
    plot.data <- data.frame(
      cell = object@meta.data$cell,
      pseudotime = object@meta.data$pseudotime,
      color.by = object@meta.data[, which(colnames(object@meta.data) == color.by)]
    )
  }

  if (!is.factor(plot.data$color.by)) plot.data$color.by <- as.factor(as.character(plot.data$color.by))

  gg <- ggplot(plot.data, aes(x=pseudotime, colour = color.by))
  color.num <- unique(plot.data$color.by)
  color.lib <- colorRampPalette(color.theme)(length(color.num))
  gg <- gg + scale_color_manual(values = color.lib)
  gg <- gg + geom_density(adjust = adjust)

  gg <- gg + theme_base()
  gg <- gg + labs(title = paste0(main))

  return(gg)
}









