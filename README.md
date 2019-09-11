
# flowSpy v1.2.7  <img src="https://github.com/JhuangLab/flowSpy/blob/master/inst/figures/logo.png" align="right" height=150 width=150/>

flowSpy is an R package to implement cellular subpopulations identification, trajectory inference, pseudotime estimation and visualization for flow and mass cytometry data. This package is developed and maintained by [JhuangLab](https://github.com/JhuangLab) at Shanghai Institute of Hematology.

Instructions and standard workflow can be found at:

https://github.com/JhuangLab/flowSpy/tree/master/vignettes

And **PDF** version of the instructions and standard workflow can be found at:

 - [**Quick_start_of_flowSpy**](https://github.com/JhuangLab/flowSpy/tree/master/inst/doc/Quick_start_of_flowSpy.pdf)

Use cases could be found at: 

https://github.com/JhuangLab/flowSpy-dataset

And **PDF** version of the specific workflows for flow and mass cytometry data can be found at:

 - [**usecase1_2**](https://github.com/JhuangLab/flowSpy-dataset/tree/master/Rmarkdown/usecase1_2.pdf) 
 - [**usecase3_4**](https://github.com/JhuangLab/flowSpy-dataset/tree/master/Rmarkdown/usecase3_4.pdf)


You can view and clone the repository of flowSpy on GitHub at:

https://github.com/JhuangLab/flowSpy


## 1 Introduction

Multidimensional single-cell-based flow and mass cytometry  enable ones to analyze multiple single-cell parameters and identify cellular populations. 
Based on classical software for analyzing [Flow Cytometry Standard](https://en.wikipedia.org/wiki/Flow_Cytometry_Standard) (FCS) data such as [`flowSOM`](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html)[1] and [`SPADE`](https://github.com/nolanlab/spade)[2], methods for inferencing cellular trajectory during a biological process are very important. 
To objectively inference differential trajectory based on time courses FCS data, we present [`flowSpy`](https://github.com/JhuangLab/flowSpy), a trajectory inference and visualization toolkit for flow and mass cytometry data. 

`flowSpy` can help you to perform four main types of analysis:

- **Clustering**. `flowSpy` can help you to discover and identify subtypes of cells. 

- **Dimensionality Reduction**. Several dimensionality reduction methods are provided in `flowSpy` package such as Principal Components Analysis (PCA), t-distributed Stochastic Neighbor Embedding (tSNE), Diffusion Maps and Uniform Manifold Approximation and Projection (UMAP). flowSpy provides both cell-based and cluster-based dimensionality reduction.

- **Trajectory Inference**. `flowSpy` can help you to construct the cellular differential based on minimum spanning tree (MST) algorithm. 

- **Pseudotime and Intermediate states definition**. The root cells need to be defined by users. The trajctroy value will be calculated based on Shortest Path from root cells and leaf cells using R `igraph` package. Subset FCS data set in `flowSpy` and find the key intermediate cell states based on trajectory value.

## 2 Installation

### 2.1 From Github

This requires the `devtools` package to be installed first.

```

# If not already installed
install.packages("devtools") 
devtools::install_github("JhuangLab/flowSpy")

library(flowSpy)

```

## 3 Workflow of flowSpy

<center> <img src="https://github.com/JhuangLab/flowSpy/blob/master/inst/figures/Workflow.png" alt="Workflow of flowSpy" /> </center>

**Workflow of flowSpy**

<center> <img src="https://github.com/JhuangLab/flowSpy/blob/master/inst/figures/algorithm.png" alt="Algorithm of flowSpy" height=70% width=70% /> </center>

**Trajectory construction and pseudotime estimation of flowSpy workflow**

## 4 Quick start (Standard Workflow)

``` {r}

# Loading packages
suppressMessages({
library(ggplot2)
library(flowCore)
library(flowSpy)
library(stringr)
})

# Read your FCS files
fcs.file <- "path to your fcs file (*.fcs)"
exp.data <- runExprsExtract(fcs.file)
# Or for more than one fcs file
fcs.data <- runExprsMerge(fcs.file)

# Build FSPY object
meta.data <- data.frame(cell = rownames(fcs.data),
                        stage = "D0" )
markers <- colnames(fcs.file)[grep("^CD", colnames(fcs.file))]
fspy <- createFSPY(raw.data = fcs.data, markers = markers, meta.data = meta.data)

# Workflow of flowSpy
fspy <- runCluster(fspy)
fspy <- processingCluster(fspy)
fspy <- runFastPCA(fspy)
fspy <- runTSNE(fspy)
fspy <- runDiffusionMap(fspy)
fspy <- runUMAP(fspy)
fspy <- buildTree(fspy, dim.type = "umap", dim.use = 1:2)
fspy <- defRootCells(fspy, root.cells = "Root cells")
fspy <- runPseudotime(fspy)
fspy <- defLeafCells(fspy, leaf.cells = "Leaf cells")
fspy <- runWalk(fspy)


```

To see the detail version, please see the vignettes [**Quick_start_of_flowSpy**](https://github.com/JhuangLab/flowSpy/tree/master/inst/doc/Quick_start_of_flowSpy.pdf)

## 5 Reported bugs and solutions

If there is any error in installing or librarying the `flowSpy` package, please contact us via e-mail forlynna@sjtu.edu.cn

## 6 Version History

Sep 11, 2019
  - Version 1.2.7
  - Changes:
    - Fixed some bugs
    - Fixed errors and warnings in `R CMD check` and `R CMD BiocCheck`

Aug 29, 2019
  - Version 1.2.6
  - Changes:
    - Fixed some bugs
    - Add trajectory plot function
    - Update vignette tutorial

Aug 14, 2019
  - Version 1.2.5
  - Changes:
    - Fixed some bugs
    - Add branch analysis and differentially expressed markers analysis

Aug 08, 2019
  - Version 1.2.4
  - Changes:
    - Fixed some bugs and finished vignette tutorial

July 24, 2019
  - Version 1.2.3
  - Changes:
    - Add phenoGraph algorithm

July 19, 2019
  - Version 1.2.2
  - Changes:
    - Fixed some bugs on cluster based downsampling


## 7 Reference

[1] Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2019). FlowSOM: Using
  self-organizing maps for visualization and interpretation of cytometry data.
  http://www.r-project.org, http://dambi.ugent.be.

[2] Qiu, P., et al., Extracting a cellular hierarchy from high-dimensional cytometry data with SPADE. Nat Biotechnol, 2011. 29(10): p.886-91.





