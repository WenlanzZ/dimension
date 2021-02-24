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
library(dplyr)
# dimension estimation
x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
results <- dimension(x, components = 1:50)

## equivelantly, if subsapce is calcualted
Subspace <- subspace(x, components = 1:50)
results  <- dimension(Subspace, method = "double_posterior")
# results <- x %>% subspace(1:50) %>% dimension()
str(results)

# truncate matrix
x_denoised <- truncate(x, components = 20, method = "threshold", alpha = 0.9, zeroout = TRUE)
x_denoised <- truncate(x, components = 20, method = "hard", zeroout = FALSE)
# equivalently, if Subspace is calculated
x_denoised <- Subspace %>% truncate(location = 1:5)
# x_denoised <- x %>% subspace(1:50) %>% truncate(location = 1:5)

# plot results
x %>% dimension() %>% plot()
x %>% dimension() %>% legacyplot(annotation =10)
```

## Using the dimension package
Click on this [link](https://rpubs.com/WenlanzZ/578132) to the vignettes for details.


## IPF single cell altas anlysis
Click on this [link](https://rpubs.com/WenlanzZ/581839) to see a workflow to include dimension in single cell RNA-Seq analysis with Seurat.
