##############################
# This is a test panel
##############################
if (F) {

  rm(list = ls())
  library(roxygen2)
  roxygenize()


verbose = T
cell.number = 500

sample.list <- paste0("D", c(0, 2, 4, 6, 8, 10))
raw <- NULL
for (i in 1:length(sample.list)) {
  sub <- read.table(paste0("../flowSpy-dataset/CleanData/usecase2/", sample.list[i], ".sub", cell.number, ".txt"), header = T, stringsAsFactors = F)
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
#save(fspy.meta.data, fspy.raw.data, file = "data/FSPYdata.rda")
batch <- factor(fspy.meta.data$stage, labels = 1:length(unique(raw$sample)))



###################################
library(ggplot2)
#library(flowSpy)
#data("FSPYdata")

markers <- c("CD34", "CD43", "CD90", "CD38", "CD49f", "CD31", "CD45RA", "FLK1", "CD73")


object <- createFSPY(raw.data = fspy.raw.data, markers = markers,
                     meta.data = fspy.meta.data,
                     normalization.method = "log",
                     verbose = T)

set.seed(1)
object <- runCluster(object, cluster.method = "som", xdim = 6, ydim = 6, verbose = T)
table(object@meta.data$cluster.id)
set.seed(2)
object <- processingCluster(object, downsampleing.size = 0.5, seed = 2)

object

object <- runKNN(object, knn = 200, verbose = T)

object <- runFastPCA(object, verbose = T)

object <- runTSNE(object, verbose = T)

object <- runDiffusionMap(object, verbose = T)

object <- runUMAP(object, verbose = T)

object <- buildTree(object, dim.type = "umap", dim.use = 1:2, verbose = T)

##########
plot2D(object, item.use = c("UMAP_1", "UMAP_2"), color.by = "stage", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)

plot2D(object, item.use = c("UMAP_1", "UMAP_2"), color.by = "som.id", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)

# pseudotime

object <- defRootCells(object, root.cells = c(13), verbose = T)

object <- runPseudotime(object, verbose = T)

object <- defLeafCells(object, leaf.cells = c(32,26,27), verbose = T)

object <- runWalk(object, verbose = T)

fetch.cells <- fetchCell(object, traj.value.log = 0.1, is.root.cells = 1, is.leaf.cells = 1)
sub.obj <- subsetFSPY(object, cells = fetch.cells)

# 4,29,9,19

plotMarkerDensity(object)


plot2D(object, item.use = c("UMAP_1", "UMAP_2"), color.by = "pseudotime", alpha = 1, main = "PCA", category = "numeric", show.cluser.id = T)

#D10.percent.stage
plotTree(object, color.by = "D6.percent", show.node.name = T, cex.size = 1.5) + scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222", "#FF3222"))

plotPieTree(object, cex.size = 1, size.by.cell.number = F) + scale_fill_manual(values = c("#00599F", "#33CC33", "#FF3222", "#F4D31D", "#FF3222","#7a06a0"))

plotPieCluster(object, item.use = c("PC_1", "PC_2"), cex.size = 0.5)
plotCluster(object, item.use = c("PC_1", "PC_2"), size = 10)


plotTree(object, color.by = "pseudotime", show.node.name = T, cex.size = 1, as.tree = T) + scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7a06a0"))



plotPseudotimeDensity(object)
plotPseudotimeTraj(object, var.cols = T) + scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7a06a0"))
plotPseudotimeTraj(object, cutoff = 0.1, var.cols = T) + scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7a06a0"))

plot2D(object, item.use = c("pseudotime", "CD43"), color.by = "som.id", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)

plot2D(object, item.use = c("CD43", "CD34"), color.by = "stage", alpha = 1, main = "PCA")

plot3D(object, item.use = c("PC1", "PC2", "PC3"), color.by = "CD43")

plot2D(object, item.use = c("tSNE_1", "tSNE_2"), color.by = "pseudotime", alpha = 1, main = "PCA") + scale_colour_gradientn(colors = c("blue", "blue", "white", "red"))

plot2D(object, item.use = c("pseudotime", "traj.value.log"), color.by = "stage")

plotTree(object, show.node.name = T, cex.size = 2) + scale_colour_gradientn(colors = "#666666")

plotTree(object, color.by = "pseudotime", as.tree = T, show.node.name = T)  + scale_colour_gradientn(colors = c("#00599F", "#00599F","#EEEEEE", "#FF3222","#FF3222"))

plot.info <- fetchPlotMeta(object, verbose = F)
ggplot(plot.info, aes(x=traj.value.log, colour = stage)) + geom_density() + theme_bw()

pdata <- aggregate(plot.info[, c(markers, "pseudotime")], list(id = plot.info$som.id), mean)
pdata <- pdata[, -1]

pdata <- pdata[order(pdata$pseudotime), ]

plot.info.sub <- plot.info[, c(markers, "pseudotime")]
plot.info.sub <- plot.info.sub[order(plot.info.sub$pseudotime), ]
plot.info.sub <- plot.info.sub[, match(c("CD34", "CD43","CD31","CD45RA","CD38","CD49f","CD90","FLK1","CD73"), colnames(plot.info.sub))]
pheatmap(t(plot.info.sub),
         color = colorRampPalette(c(rep("#00599F",3), "#FFFFFF", rep("#FF3222",3)))(100),
         cluster_rows = F, cluster_cols = F, scale = "row", fontsize_col = 0.01)


plot.info.sub <- plot.info[, c("pseudotime", "traj.value.log")]
plot.info.sub <- plot.info.sub[plot.info.sub$traj.value.log > 0.2, ]
plot.info.sub <- plot.info.sub[order(plot.info.sub$pseudotime),]
plot(plot.info.sub$pseudotime, ylim = c(0,1), xlim = c(0,12000), col = "#AA00C2")

###############
sub.obj <- runKNN(sub.obj, knn = 30)

set.seed(1)
sub.obj <- runCluster(sub.obj, cluster.method = "som", xdim = 2, ydim = 3)
table(sub.obj@meta.data$cluster.id)

sub.obj <- runFastPCA(sub.obj)
sub.obj <- runTSNE(sub.obj)
sub.obj <- runDiffusionMap(sub.obj)
sub.obj <- runUMAP(sub.obj)
sub.obj <- buildTree(sub.obj, cluster.type = "som", dim.type = "tSNE")

plot(sub.obj@network$mst)

plot2D(sub.obj, item.use = c("CD34", "CD49f"), color.by = "stage", alpha = 1, main = "PCA", category = "categorical", show.cluser.id = T)

plotTree(sub.obj, color.by = "D8.percent", show.node.name = T, cex.size = 2) + scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))

plotPseudotimeDensity(sub.obj, adjust = 1)

run.time <- read.xlsx("../dataset/RunningTime.xlsx")
p <- ggbarplot(run.time, x = "Function", y = "Time_5000", fill = "#0076D4",
               xlab = "", ylab = "Time (Seconds)", label = T) + theme_base()
p <- p + scale_y_continuous(limits = c(0,270), expand=c(0,0)) + coord_flip()
p <- p + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))
ggsave("../dataset/Time_5000.pdf", width = 6, height = 5)


###############
# merge FCS
fcsFiles <- paste0("../flowSpy-dataset/FCS/usecase2/D", c(0,2,4,6,8,10), ".fcs")

merged <- runExprsMerge(fcsFiles)



}


















