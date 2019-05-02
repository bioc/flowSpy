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
  sub <- read.table(paste0("inst/extdata/dataset/", sample.list[i], ".sub.txt"), header = T, stringsAsFactors = F)
  sub$sample <- sample.list[i]
  sub <- sub[1:cell.number, ]
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



object <- createFSPY(raw.data = raw.data, markers = markers,
                     meta.data = meta.data,
                     log.transformed = F,
                     verbose = T)

object <- runCluster(object, cluster.method = "som", xdim = 6, ydim = 6)

object <- runFastPCA(object)

object <- runTSNE(object)

object <- runDiffusionMap(object)

object <- runUMAP(object)

object <- runSOM(object)

object <- runMclust(object)


xdim = 4
ydim = 4
rlen = 8
mst = 1
alpha = c(0.05, 0.01)
radius = 1
init = FALSE
distf = 2
codes = NULL
importance = NULL
method = "euclidean"
verbose= T

object <- runSOM(object, xdim = xdim, ydim = ydim)

object <- updatePlotMeta(object)

item.use = c("tSNE1", "tSNE2")
color.by = "stage"
order.by = NULL
size = 1
alpha = 1
category = NULL
main = "2D plot of FSPY"
plot.theme = theme_base()
trajectory = "som"
show.node.name = T
color.theme = NULL

############ test of plot
plot2D(object, item.use = c("PC1", "PC2"), color.by = "mc.id", alpha = 0.6, main = "PCA", category = "categorical")
plot2D(object, item.use = c("tSNE1", "tSNE2"), color.by = "mc.id", alpha = 0.6, main = "tSNE", category = "categorical", trajectory = "")
plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "som.id", alpha = 1, main = "UMAP", category = "categorical", trajectory = "som")

plot2D(object, item.use = c("UMAP1", "UMAP2"), color.by = "pseudotime", alpha = 0.6, main = "PCA")
plot3D(object, item.use = c("PC1", "PC2", "PC3"), color.by = "CD34", size = 0.5,
       angle = 45, main = "pseudotime")
plotSOM(object, color.by = "CD34", show.node.name = T)
plotSOMtree(object, color.by = "trunk.id",
            show.node.name = T, cex.size = 1.5,
            color.theme = "#FFCC66")


data("slingshotExample")
sds <- slingshot(rd, clusterLabels = cl, start.clus = '1')


id <- as.vector(object@meta.data$som.id)
names(id) <- rownames(object@umap.layout[, 1:2])
sds <- slingshot(object@pca.load[, 1:6], clusterLabels = object@meta.data$som.id, start.clus = '11')

plot(object@tsne.value[, 1:2], col = object@log.data[, 2], pch = 16)
lines(sds, lwd = 3)

return(object)





############ DDRtree
library(DDRTree)

a <- DDRTree(t(object@umap.layout[, 1:2]), dimensions = 2)

Z <- a$Z
Y <- a$Y
W <- a$W
stree <- a$stree

plot(Z[1, ], Z[2, ])
plot(Y[1, ], Y[2, ])
plot(W[, 1], W[, 2])


fullGraph <- igraph::graph.adjacency(as.matrix(stree),
                                     mode = "undirected",
                                     weighted = TRUE)
plot(fullGraph)

save(object, file = "0328.fspy.Robj")
load("0328.fspy.Robj")


mat <- t(object@log.data)
adj <- matrix(0, ncol(mat), ncol(mat))
rownames(adj) <- colnames(adj) <- colnames(mat)
for(i in seq_len(ncol(mat))) {
  adj[i, colnames(mat)[object@knn.index[i,]]] <- 1
}
g <- igraph::graph.adjacency(adj, mode="undirected")
# remove self loops
g <- simplify(g)
## identify communities
km <- igraph::cluster_fast_greedy(g)
km <- igraph::cluster_spinglass(g)
sizes(km)

root.cells <- which(object@meta.data$branch.id == "3-3")
walk <- random_walk(object@network$knn.G, start = root.cells, step = 12000)
walk



#ggsave("1.plotGATE.dotmesh.pdf", p, width = 6, height = 5)

########## KNN cluster
library(RANN)

mat <- t(object@gate.data)

knn.info <- RANN::nn2(object@gate.data, k=25)

knn <- knn.info$nn.idx
adj <- matrix(0, ncol(mat), ncol(mat))
rownames(adj) <- colnames(adj) <- colnames(mat)
for(i in seq_len(ncol(mat))) {
  adj[i,colnames(mat)[knn[i,]]] <- 1
}

library(igraph)
g <- igraph::graph.adjacency(adj, mode="undirected")
g <- simplify(g)
## identify communities
km <- igraph::cluster_walktrap(g)
com <- km$membership
names(com) <- km$names

## compare to annotations
object@meta.data$knn.id <- com


object <- runFastPCA(object)

object <- runTSNE(object)

# fast tsne

source('../FIt-SNE-master/fast_tsne.R', chdir=T)
fltsne <- fftRtsne(object@gate.data, fast_tsne_path = "/Users/daiyuting/Documents/projects/bioVis/flowSpy/v1/FIt-SNE-master/bin/fast_tsne")


object <- runDiffusionMap(object)

a <- umap(object@gate.data, config = umap.defaults)

object@pca.load <- cbind(object@pca.load, a$layout)
colnames(object@pca.load) <- c(colnames(object@pca.load)[1:9], "PC01", "PC02")

object@pca.load[, 10] <- a$layout[, 1]
object@pca.load[, 11] <- a$layout[, 2]

kmeams.cluster <- kmeans(object@gate.data, 25, nstart = 20)

table(kmeams.cluster$cluster)

object@meta.data$som.node.id <- kmeams.cluster$cluster

plot2D(object, item.use = c("PC1", "PC2"), color.by = "stage", alpha = 0.6, main = "PCA", size = 0.8)
plot2D(object, item.use = c("PC01", "PC02"), color.by = "knn.id", alpha = 0.3, main = "PCA", size = 0.8)
plot2D(object, item.use = c("PC01", "PC02"), color.by = "stage", alpha = 0.3, main = "PCA", size = 0.8)
plot2D(object, item.use = c("tSNE_1", "tSNE_2"), color.by = "knn.id", alpha = 1, main = "tSNE", size = 0.8)
plot2D(object, item.use = c("DC1", "DC2"), color.by = "stage", alpha = 0.6, main = "Diffusion Map", size = 0.8)
plot3D(object, item.use = c("DC1", "DC2","DC3"), color.by = "stage", size = 0.5, angle = 45, main = "Diffusion Map")
plot3D(object, item.use = c("PC1", "PC2","PC3"), color.by = "stage", size = 0.5, angle = 45, main = "PCA")


# paramter of som
xdim = 10
ydim = 10
rlen = 10
mst = 1
alpha = c(0.05,  0.01)
radius = 1
init = FALSE
distf = 2
silent = FALSE
codes = NULL
importance = NULL
method = "euclidean"
verbose= T

object <- runSOM(object, xdim = 5, ydim = 5)

###### plot SOM
plotSOM(object, color.by = "cell.percent",
        color.theme = c("#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))
plotSOM(object, color.by = "CD49f", show.node.name = T, cex.size = 1.2,
        color.theme = c("#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))

plotSOM(object, color.by = "D10.som.percent", show.node.name = T, cex.size = 1.2,
        color.theme = c("#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))
plotSOM(object, color.by = "D4.som.percent", show.node.name = T, cex.size = 1.2,
        color.theme = c("#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))
plotSOM(object, color.by = "D0.som.percent", show.node.name = T, cex.size = 1.2,
        color.theme = c("#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))

# pseudotime
root.cell <- as.character(object@meta.data$cell[which(object@meta.data$som.node.id == 2)])

root.cell <- sample(root.cell, size = floor(length(root.cell)/4))
root.cell <- root.cell[grep("D0", root.cell)]

object <- defRootCells(object, root.cell = root.cell)
object <- pseudotimeProcess(object)

plot2D(object, item.use = c("DC1", "DC2"), color.by = "knn.id", alpha = 0.9,
       main = "pseudotime", size = 0.5,
       color.theme = c("#f2de00", "#8bd129", "#33CCFF", "#0066CC", "#CC66FF", "#FF3300"))

plot2D(object, item.use = c("PC01", "PC02"), color.by = "CD43", alpha = 0.9, main = "pseudotime", size = 0.5, color.theme = c("#0000CC", "#0000CC", "#0000CC", "#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))

plot3D(object, item.use = c("PC01", "PC02", "pseudotime"), color.by = "knn.id", size = 0.5,
       angle = 45, main = "pseudotime", color.theme = c("#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))

plot2D(object, item.use = c("tSNE_1", "tSNE_2"), color.by = "knn.id", alpha = 0.9, main = "stage", size = 0.5)

plot3D(object, item.use = c("pseudotime", "DC1", "DC2"), color.by = "stage", size = 0.5,
       angle = 60, main = "stage")


plotPseudotimeDensity(object, color.by = "knn.id")


#D10.som.percent
plotSOMtree(object, color.by = "aa",
            show.node.name = T, cex.size = 1.5,
            color.theme = "#FFCC66")

plotSOMtree(object, color.by = "pseudotime",
            show.node.name = T, cex.size = 1.5,
            color.theme = c("#0000CC", "#0000CC", "#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))

plotSOMtree(object, color.by = "CD90", show.node.name = T, cex.size = 1.5)

object@meta.data$som.target = 0
som.id =  c(2,16)
for (i in som.id) {
  object@meta.data$som.target[which(object@meta.data$som.node.id == i)] = i
}
plotGATE(object, plot.markers = c("CD43", "CD49f"), color.by = "som.target",
         plot.type = "dot", alpha = 0.5,
         color.theme = c("#666666", rainbow(length(som.id))))
plotPseudotimeDensity(object, color.by = "som.target")

object@meta.data$som.target <- as.factor(object@meta.data$som.target)

plot2D(object, item.use = c("tSNE_1", "tSNE_2"),
       color.by = "som.target", alpha = 0.6, main = "tSNE", size = 0.5,
       color.theme = colorRampPalette(c("#CCCCCC", "#0000CC", "#0099FF", "#339900", "#FFFF00", "#FF6600", "#FF0000"))(11)  )


som.id =  c(31)
a <- object@meta.data[ object@meta.data$som.target %in% som.id, ]
aa <- object@gate.data[ rownames(object@gate.data) %in% a$cell, ]
summary(aa)


plot.sub <- data.frame(
  object@meta.data,
  raw.data[match(object@meta.data$cell, rownames(raw.data)), ]
)
ggscatter(plot.sub, x = "pseudotime", y = "CD34", size = 1)


som.id <- c(16, 2)
mat <- object@som$node.attr
som.sub.mat <- mat[rownames(mat) %in% som.id, ]
id.mat <- object@meta.data[object@meta.data$som.node.id %in% id.som, ]

plot.mat <- matrix(0, nrow = length(id.som), ncol = ncol(raw.data))
rownames(plot.mat) <- id.som
colnames(plot.mat) <- colnames(raw.data)

for (i in 1:nrow(plot.mat)) {
  sub.1 <- as.character(id.mat$cell[which(id.mat$som.node.id == rownames(plot.mat)[i])])
  plot.mat[i, ] <- colMedians(raw.data[match(sub.1, rownames(raw.data)), ])
}
markers.plot <- c("Brachury", "CD45RA", "CD43", "OCT4_", "CD38", "CD90", "CD31", "CD49f", "CD73", "CD45", "FLK1", "CD34")
plot.mat.this <- t(plot.mat[, markers.plot])
plot.mat.this[2, ] = plot.mat.this[2, ] - 1
plot.mat.this[1, ] = plot.mat.this[1, ] - 0.5
#plot.mat.this <- plot.mat.this - rowMedians(plot.mat.this)
pheatmap(plot.mat.this, cluster_cols = F,
         cluster_rows = T, scale = "none",
         border_color = NA,
         color = colorRampPalette( c("#3399FF", "#FFFFFF", "#FFCC33", "#FF6600") )(100))


library(roxygen2)
roxygenize()


############# other plot
idx <- sample(1:dim(object@dm@transitions)[1], size = 500)
idx <- 1:500
trans.mat <- as.matrix(object@dm@transitions[idx, idx])

p <- pheatmap(trans.mat, fontsize = 0.1)
ggsave("3.trans.mat.pdf", p, width = 10, height = 9)






library(dbscan)

data(iris)
iris <- as.matrix(iris[,1:4])
kNNdistplot(iris, k = 5)
abline(h=.5, col = "red", lty=2)

res <- dbscan(iris, eps = .5, minPts = 5)
res
pairs(iris, col = res$cluster + 1L)


}


















