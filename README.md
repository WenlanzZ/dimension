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
x <- x_sim(n = 100, p = 150, ncc = 30, var = c(rep(10,5),rep(2,25)))
results <- dimension(x, components = 1:50, times = 10, p = 0.95)

## equivelantly, if subsapce is calcualted
Subspace <- subspace(x, components = 1:50, times = 10)
results  <- dimension(subspace_ = Subspace)
str(results)

# clip matrix
x_clp <- clipped(x, components = 20, method = "threshold", alpha = 0.9, zeroout = TRUE)
x_clp <- clipped(x, components = 20, method = "hard", zeroout = FALSE)
# equivalently, if Subspace is calculated
x_clp <- clipped(subspace_ = Subspace, method = "identity", location = c(1:5))
x_clp

# plot results
plot(Subspace, annotation = 30, changepoint = results$dimension)
modified_legacyplot(results$Changepoint$bcp_irl, annotation = 30)
modified_legacyplot(results$Changepoint$bcp_post, annotation = 30)

## IPF single cell altas anlysis
Click on this [link](dimension.html) to see a workflow to include dimension in single cell RNA-Seq analysis with Seraut.