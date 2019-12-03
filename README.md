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
## basic example code
X <- Xsim(n = 150, p = 100, ncc = 10, var = 2)
results <- dimension(X, rank = 1:40, times = 10, basis="eigen")
 
#equivelantly, if subsapce is calcualted
Subspace <- subspace(X, rank = 1:40, times = 10,  basis = "eigen")
results <- dimension(subspace_ = Subspace)
str(results)
plot(results$Subspace, Changepoint = results$Changepoint$dimension, annotation = 10)
modified_legacyplot(results$Changepoint$bcp_irl, annotation = 10)
modified_legacyplot(results$Changepoint$bcp_post, annotation = 10)
```
