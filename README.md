# dimension

<!-- badges: start -->
<!-- badges: end -->

The goal of "dimension" is to estimate the dimension of a signal-rich subspace in large high-dimensional data.

## Installation

You can install the released version of dimension from GitHub with:

``` r
devtools::install_github("WenlanzZ/dimension")
```

## Example

This is a basic example which shows you how to use this package:

``` r
library(dimension)
# dimension estimation
x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
results <- dimension(x, components = 1:50)

## equivelantly, if subsapce is calcualted
Subspace <- subspace(x, components = 1:50)
results  <- dimension(subspace_ = Subspace)
str(results)

# truncate matrix
x_denoised <- truncate(x, components = 20, method = "threshold", alpha = 0.9, zeroout = TRUE)
x_denoised <- truncate(x, components = 20, method = "hard", zeroout = FALSE)
# equivalently, if Subspace is calculated
x_denoised <- truncate(subspace_ = Subspace, method = "identity", location = c(1:5))
x_denoised

# plot results
plot(results$Subspace, changepoint = results$dimension, annotation = 10)
modified_legacyplot(results$bcp_irl, annotation = 10)
```

## Using the dimension package
Click on this [link](https://rpubs.com/WenlanzZ/578132) to the vignettes for details.


## IPF single cell altas anlysis
Click on this [link](https://rpubs.com/WenlanzZ/581839) to see a workflow to include dimension in single cell RNA-Seq analysis with Seurat.
