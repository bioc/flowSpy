# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


if (F) {

rm(list = ls())
library(roxygen2)
roxygenize()


verbose = T
cell.number = 500

sample.list <- paste0("D", c(0, 2, 4, 6, 8, 10))
raw <- NULL
for (i in 1:length(sample.list)) {
  sub <- read.table(paste0("../dataset/", sample.list[i], ".sub.txt"), header = T, stringsAsFactors = F)
  sub$sample <- sample.list[i]
  #sub <- sub[sample(1:dim(sub)[1], cell.number), ]
  raw <- rbind(raw, sub)
}
table(raw$sample)
raw$sample <- factor(as.character(raw$sample), levels = sample.list)



markers <- c("CD34", "CD43", "CD38", "CD90", "CD49f", "CD31", "CD45RA", "FLK1", "CD73")
length(markers)

raw.data <- as.matrix(raw[, 1:(dim(raw)[2]-1) ])
rownames(raw.data) <- paste0(raw$sample, "_", 1:length(raw$sample))
meta.data <- data.frame(cell = paste0(raw$sample, "_", 1:length(raw$sample)),
                        stage = raw$sample)

fspy.meta.data <- meta.data
fspy.raw.data <- raw.data

#save(fspy.meta.data, fspy.raw.data, file = "data/fspy.data.rda")

object <- createFSPY(raw.data = fspy.raw.data, markers = markers,
                     meta.data = fspy.meta.data,
                     log.transform = F,
                     verbose = T)


object <- runKNN(object, knn = 30)

set.seed(1)
object <- runCluster(object, cluster.method = "som", xdim = 6, ydim = 6)
table(object@meta.data$cluster.id)


object <- runFastPCA(object)

object <- runTSNE(object)

object <- runDiffusionMap(object)

object <- runUMAP(object)

object <- updatePlotMeta(object)

object <- buildTree(object, cluster.type = "som", dim.type = "umap")

plot(object@network$mst)

# pseudotime
object <- defRootCells(object, root.cells = c(24))

object <- runPseudotime(object)

object <- defLeafCells(object, leaf.cells = c(6,29,32,35), pseudotime.cutoff = 0.5)

object <- runWalk(object)


plotPseudotimeDensity(object)
plotPseudotimeTraj(object, var.cols = T) + scale_colour_gradientn(colors = c("blue", "red"))
plotPseudotimeTraj(object, cutoff = 0.2, var.cols = T) + scale_colour_gradientn(colors = c("blue", "red"))





plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "som.id", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)

plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "stage", alpha = 1, main = "PCA", category = "categorical")

plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "traj.value.log", alpha = 0.5, main = "PCA", category = "numeric") + scale_colour_gradientn(colors = c("#FFFFCC", "red", "red", "red"))

plot3D(object, item.use = c("UMAP1", "UMAP2", "pseudotime"), color.by = "stage")

plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "CD49f", alpha = 1, main = "PCA") + scale_colour_gradientn(colors = c("blue", "blue", "white", "red"))

plot2D(object, item.use = c("pseudotime", "traj.value.log"), color.by = "stage")

plotTree(object, color.by = "CD49f", as.tree = T, show.node.name = T, root.id = 24)  + scale_colour_gradientn(colors = c("blue", "white", "red"))

plot.info <- fetchPlotMeta(object, verbose = F)
ggplot(plot.info, aes(x=pseudotime, colour = stage)) + geom_density() + theme_base()




save(object, file = "0505.FSPY.Robj")




}


















