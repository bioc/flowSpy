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
cell.number = 2000

sample.list <- paste0("D", c(0, 2, 4, 6, 8, 10))
raw <- NULL
for (i in 1:length(sample.list)) {
  sub <- read.table(paste0("../dataset/", sample.list[i], ".sub", cell.number, ".txt"), header = T, stringsAsFactors = F)
  sub$sample <- sample.list[i]
  raw <- rbind(raw, sub)
}
raw$sample <- factor(as.character(raw$sample), levels = sample.list)
table(raw$sample)

markers <- c("CD34", "CD43", "CD38", "CD90", "CD49f", "CD31", "CD45RA", "FLK1", "CD73")
length(markers)

raw.data <- as.matrix(raw[, 1:(dim(raw)[2]-1) ])
rownames(raw.data) <- paste0(raw$sample, "_", 1:length(raw$sample))
meta.data <- data.frame(cell = paste0(raw$sample, "_", 1:length(raw$sample)),
                        stage = raw$sample)

fspy.meta.data <- meta.data
fspy.raw.data <- raw.data

save(fspy.meta.data, fspy.raw.data, file = "data/FSPYdata.rda")


batch <- factor(fspy.meta.data$stage, labels = 1:length(unique(raw$sample)))

object <- createFSPY(raw.data = fspy.raw.data, markers = markers,
                     meta.data = fspy.meta.data,
                     log.transform = T,
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

plot2D(object, item.use = c("tSNE1", "tSNE2"), color.by = "som.id", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)

plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "stage", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)

plot2D(object, item.use = c("DC1", "DC2"), color.by = "stage", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)



# pseudotime
object <- defRootCells(object, root.cells = c(18,11))

object <- runPseudotime(object)

object <- defLeafCells(object, leaf.cells = c(5,12,29), pseudotime.cutoff = 0.5)

object <- runWalk(object)


plotPseudotimeDensity(object)
plotPseudotimeTraj(object, var.cols = T) + scale_colour_gradientn(colors = c("#00599F",  "#EEEEEE", "#FF3222"))
plotPseudotimeTraj(object, cutoff = 0.4, var.cols = T) + scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))


plot2D(object, item.use = c("pseudotime", "CD43"), color.by = "som.id", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)



plot2D(object, item.use = c("CD43", "CD34"), color.by = "stage", alpha = 1, main = "PCA")

plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "traj.value.log", alpha = 0.5, main = "PCA", category = "numeric") + scale_colour_gradientn(colors = c("#FFFFCC", "red", "red", "red"))

plot3D(object, item.use = c("UMAP1", "UMAP2", "pseudotime"), color.by = "stage")

plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "CD49f", alpha = 1, main = "PCA") + scale_colour_gradientn(colors = c("blue", "blue", "white", "red"))

plot2D(object, item.use = c("pseudotime", "traj.value.log"), color.by = "stage")

plotTree(object, show.node.name = T, cex.size = 2) + scale_colour_gradientn(colors = "#666666")

plotTree(object, color.by = "pseudotime", as.tree = T, show.node.name = T, root.id = 24)  + scale_colour_gradientn(colors = c("#00599F", "#00599F","#EEEEEE", "#FF3222","#FF3222"))

plot.info <- fetchPlotMeta(object, verbose = F)
ggplot(plot.info, aes(x=traj.value.log, colour = stage)) + geom_density() + theme_bw()

pdata <- aggregate(plot.info[, c(markers, "pseudotime")], list(id = plot.info$som.id), mean)
pdata <- pdata[, -1]

pdata <- pdata[order(pdata$pseudotime), ]

plot.info.sub <- plot.info[, c(markers, "pseudotime")]
plot.info.sub <- plot.info.sub[order(plot.info.sub$pseudotime), ]
plot.info.sub <- plot.info.sub[, match(c("CD34", "CD43","CD31","CD45RA","CD38","CD49f","CD90","FLK1","CD73"), colnames(plot.info.sub))]
pheatmap(plot.info.sub,
         color = colorRampPalette(c("#00599F", "#00599F", "#FFFFFF", "#FF3222", "#FF3222"))(100),
         cluster_rows = F, cluster_cols = F, scale = "column")


plot.info.sub <- plot.info[, c("pseudotime", "traj.value.log")]
plot.info.sub <- plot.info.sub[plot.info.sub$traj.value.log >= 0, ]
plot.info.sub <- plot.info.sub[order(plot.info.sub$pseudotime),]
plot(plot.info.sub$pseudotime, ylim = c(0,1))


}


















