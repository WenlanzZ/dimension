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
X <- Xsim(n = 150, p = 100, ncc = 10, var = 2)
results <- dimension(X, rank = 1:40, times = 10, basis="eigen")

## equivelantly, if subsapce is calcualted
Subspace <- subspace(X, rank = 1:40, times = 10,  basis = "eigen")
results <- dimension(subspace_ = Subspace)

# clip matrix
X_clp <- clipped(X, rank = 20, method = "threshold", alpha = 0.9, zeroout = TRUE)
x_clp<-clipped(x, rnk = 20, method = "hard", zeroout = FALSE)
x_clp<-clipped(x, rnk = 20, method = "identity", location = c(1:15), zeroout = FALSE)

## equivelantly, if Subspace is calcualted
X_clp <- clipped(subspace_ = Subspace, method = "identity", location = c(1:5), zeroout = TRUE)

# plot results
str(results)
plot(results$Subspace, Changepoint = results$Changepoint$dimension, annotation = 10)
modified_legacyplot(results$Changepoint$bcp_irl, annotation = 10)
modified_legacyplot(results$Changepoint$bcp_post, annotation = 10)
```
