
if (F) {

library(flowCore)
library(LSD)
library(stringr)


data.dir <- "/Users/daiyuting/Documents/projects/bioVis/flowSpy/data/2019-02-13/"

data.name <- "20190123_D8 VEGF S_022.fcs"
out.name = "D8_VEGF_S"

data <- flowCore::read.FCS(filename = paste0(data.dir, data.name))

recol.data <- data@parameters@data
recol.data$desc[which(is.na(recol.data$desc))] <- str_replace_all(recol.data$name, "-", "_")[which(is.na(recol.data$desc))]

##### compensation
comp.mat <- description(data)$SPILL
data.comp <- compensate(data, comp.mat)

fcs.data <- as.data.frame(exprs(data.comp))
colnames(fcs.data) <- recol.data$desc

fcs.data <- fcs.data[which((fcs.data$FSC_A > 0) & (fcs.data$FSC_A > 0) ), ]
fcs.data <- fcs.data[which((fcs.data$FSC_A > 50000) & (fcs.data$FSC_A < 200000) ), ]
fcs.data <- fcs.data[which((fcs.data$SSC_A > 30000) & (fcs.data$SSC_A < 150000) ), ]
fcs.data <- fcs.data[which((fcs.data$FSC_W > 80000) & (fcs.data$FSC_W < 130000) ), ]
fcs.data <- fcs.data[which((fcs.data$FSC_H > 50000) & (fcs.data$FSC_H < 120000) ), ]
fcs.data <- fcs.data[which((fcs.data$SSC_W > 10000) & (fcs.data$SSC_W < 100000) ), ]
fcs.data <- fcs.data[which((fcs.data$SSC_H > 10000) & (fcs.data$SSC_H < 90000) ), ]

colnames(fcs.data)

markers.1 <- "CD43"
markers.2 <- "CD34"

a <- fcs.data[, which(colnames(fcs.data) == markers.1)]
b <- fcs.data[, which(colnames(fcs.data) == markers.2)]
heatscatter(a, b, cexplot = 0.2,
            main = paste0("Flow Cytometry"), xlab = markers.1, ylab = markers.2)

colnames <- c("FSC_A", "FSC_H", "FSC_W", "SSC_A", "SSC_H", "SSC_W",
              "CD34", "CD43", "CD38", "CD90", "CD49f", "CD31", "CD45RA", "FLK1", "CD73")
colnames[!colnames %in% colnames(fcs.data)]


############## output data
colnames(fcs.data)[which(colnames(fcs.data) == "PE_Cy7_A")] = "CD38"
colnames[!colnames %in% colnames(fcs.data)]

fcs.data <- fcs.data[, match(colnames, colnames(fcs.data))]
dim(fcs.data)

write.table(fcs.data, paste0(data.dir, "dealt/", out.name, ".txt"), sep = "\t", quote = F, row.names = F)


colnames(fcs.data)[which(colnames(fcs.data) == "APC_A")] = "CD34"
colnames(fcs.data)[which(colnames(fcs.data) == "FITC_A")] = "CD43"
colnames(fcs.data)[which(colnames(fcs.data) == "PE_Cy7_A")] = "CD38"
colnames(fcs.data)[which(colnames(fcs.data) == "BV421_A")] = "CD90"
colnames(fcs.data)[which(colnames(fcs.data) == "BV650_A")] = "CD49f"
colnames(fcs.data)[which(colnames(fcs.data) == "BV605_A")] = "CD31"
colnames(fcs.data)[which(colnames(fcs.data) == "BV510_A")] = "CD45RA"
colnames(fcs.data)[which(colnames(fcs.data) == "PE_A")] = "FLK1"
colnames(fcs.data)[which(colnames(fcs.data) == "BV 735_A")] = "CD73"

# for D0
fcs.data <- fcs.data[which(fcs.data$CD43 < 10000), ]
fcs.data <- fcs.data[which(fcs.data$CD34 < 10000), ]


#######################################
#### Log transformed of FCS data set

name="D12"
fcs.data <- read.table( paste0(data.dir, "dealt/", name, ".txt"), sep = "\t", header = T )
dim(fcs.data)

markers <- c("CD34", "CD43", "CD38", "CD90", "CD49f", "CD31", "CD45RA", "FLK1", "CD73")
markers.idx <- match(markers, colnames(fcs.data))
fcs.data.log <- fcs.data

fcs.data.log[, markers.idx] <- log10(abs(fcs.data.log[, markers.idx]) + 1)

markers.1 <- "CD43"
markers.2 <- "CD49f"

a <- fcs.data.log[, which(colnames(fcs.data.log) == markers.1)]
b <- fcs.data.log[, which(colnames(fcs.data.log) == markers.2)]
heatscatter(a, b, cexplot = 0.2,
            main = paste0("Flow Cytometry"), xlab = markers.1, ylab = markers.2)

write.table(fcs.data.log, paste0(data.dir, "dealt/", name, ".log.txt"), sep = "\t", quote = F, row.names = F)


cell.size = 2000
sample.list <- c("D0", "D2", "D4", "D6", "D8", "D10", "D12")

for (i in 1:length(sample.list)) {
  sub <- read.table(paste0(data.dir, "dealt/", sample.list[i], ".log.txt"),  sep = "\t", header = T)
  idx <- sample(1:dim(sub)[1], cell.size)
  sub <- sub[idx, ]
  write.table(sub, paste0("/Users/daiyuting/Documents/projects/bioVis/flowSpy/v1/flowSpy/inst/extdata/dataset/", sample.list[i], ".sub.txt"), sep = "\t", quote = F, row.names = F)
}

}
















