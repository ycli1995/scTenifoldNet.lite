
# scTenifoldNet.lite

[![License: GPL (>=2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D2%29-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

The goal of `scTenifoldNet.lite` is to provide fast and memory-efficient 
implementation of the [`scTenifoldNet`](https://github.com/cailab-tamu/scTenifoldNet) 
workflow for constructing and comparing single-cell gene regulatory networks. 
The computational performance for principal component regression, tensor decomposition, 
and manifold alignment was optimized, while preserving the core analytical framework 
of the original `scTenifoldNet` method. 

## Installation

You can install the development version of scTenifoldNet.lite like so:

``` r
devtools::install_github("ycli1995/scTenifoldNet.lite")
```

## Package Overview

The original six-step pipeline in [`scTenifoldNet`](https://github.com/cailab-tamu/scTenifoldNet) 
has been optimized for memory and speed. Below describes the main differences.

### `scQC`: Quality control
This step filters cells by library size, outlier detection, minimum gene 
expression fraction, and mitochondrial read ratio. 

* I add `drop = FALSE` for subsetting to keep the matrix formats.

### `cpmNormalization`: Counts-per-million (CPM) normalization
* Use `RcppEigen` to speed up.

### `makeNetworks`: Constructs gene regulatory networks from subsampled cells
* The GRNs are constructed using principal component regression (`pcNet`).
* I add a `fast = TRUE` parameter to use global principal component scores among
all genes, which avoids running `RSpectra::svds()` upon each gene to speed up.

### `tensorDecomposition`: Network denoising
* I only keep the 3D CANDECOMP/PARAFAC (CP) tensor decomposition path.
* The denoised tensor will be represented as `U_list` instead of a memeory-expensive 
array.

### `manifoldAlignment`: Non-linear manifold alignment
* The `W` matrix with `2n * 2n` elements will not be constructed to save memeory.

### `dRegulation`: Differential regulation testing
* This step testing the differential regulation on the manifold alignment output.
* No need to optimize.

## Example
The motivation of this package is to perform `scTenifoldKnk()` pipeline for 
virtually knockout:

``` r
library(scTenifoldNet.lite)

scTenifoldKnk(countMatrix = scRNAseq, gKO = "G100", qc_minLSize = 0, fast = TRUE)
```

## Citation

If you use `scTenifoldNet.lite` in your research, please cite the originial manuscript:

> Osorio, D., Zhong, Y., Li, G., Huang, J. Z., & Cai, J. J. (2020). scTenifoldNet: A Machine Learning Workflow for Constructing and Comparing Transcriptome-wide Gene Regulatory Networks from Single-Cell Data. *Patterns*, **1**(9), 100139. [doi:10.1016/j.patter.2020.100139](https://doi.org/10.1016/j.patter.2020.100139)

> Osorio, Daniel, et al. "scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation." *Patterns* 3.3 (2022).
