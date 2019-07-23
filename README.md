
# flowSpy v1.2.2  <img src="https://github.com/ytdai/flowSpy/blob/master/inst/figures/logo.png" align="right" height=10 width=10/>

flowSpy is a trajectory inference and visualization toolkit of flow and mass cytometry data. This package is developed and maintained by [JhuangLab](https://github.com/JhuangLab) at Shanghai Institute of Hematology.

Instructions, documentation, and tutorials can be found at:

https://github.com/ytdai/flowSpy/tree/master/vignettes

You can view and clone the repository of flowSpy on GitHub at:

https://github.com/ytdai/flowSpy

## 1 Introduction

High throughput cell-based assays with flow cytometric signals enable ones to analyze multiple single-cell parameters and identify cellular populations. 
Based on classical software for analyzing [Flow Cytometry Standard](https://en.wikipedia.org/wiki/Flow_Cytometry_Standard) (FCS) data such as [`flowSOM`](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html)[1] and [`SPADE`](https://github.com/nolanlab/spade)[2], methods for inferencing cellular trajectory during a biological process are very important. 
To objectively inference differential trajectory based on time courses FCS data, we present [`flowSpy`](https://github.com/ytdai/flowSpy), a trajectory inference and visualization toolkit of FCS data. In this tutorial, we will present how to complete an analysis workflow of time courses FCS data using `flowSpy` package. 

`flowSpy` can help you to perform four main types of analysis:

- **Clustering**. `flowSpy` can help you to discover and identify subtypes of cells. 

- **Reducing Dimensions**. Several dimensionality reduction methods are provided in `flowSpy` package such as Principal Components Analysis (PCA), t-distributed Stochastic Neighbor Embedding (tSNE), Diffusion Maps and Uniform Manifold Approximation and Projection (UMAP).

- **Trajectory Inference**. `flowSpy` can help you to construct the cellular differential based on Minimum Spanning Tree (MST) algorithm. 

- **Pseudotime and Intermediate states definition**. The root cells need to be defined by users. The trajctroy value will be calculated based on Shortest Path from root cells and leaf cells using R `igraph` package. Subset FCS data set in `flowSpy` and find the key intermediate cell states based on trajectory value.

## 2 Installation

`flowSpy` can be installed in one of two ways:

### 2.1 From Bioconductor 

The `flowSpy` package is not uploaded to Bioconductor server yet, so now you can only download it through Github.

`flowSpy` runs in the [R statistical computing environment](https://www.r-project.org/). You will need R version 3.5 or higher to have access to the latest features. 

```

# The flowSpy package is not uploaded to Bioconductor server yet
# Please install it through GitHub.
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("flowSpy")

library(flowSpy)

```
### 2.2 From Github

This requires the `devtools` package to be installed first.

```

# If not already installed
install.packages("devtools") 
devtools::install_github("ytdai/flowSpy")

library(flowSpy)

```

## 3 Workflow of flowSpy

<center> <img src="https://github.com/ytdai/flowSpy/blob/master/inst/figures/Workflow.png" alt="Workflow of flowSpy" /> </center>

## 4 Version History

July 19, 2019
  - Version 1.2.2
  - Changes:
    - Fixed some bugs on cluster based downsampling





