# Residuals and Whitened Residuals of a Multivariate Gaussian

Computes the centered response \\r_i = y_i - \mu\\ and its image under
the inverse covariance, \\w_i = \Sigma^{-1} r_i\\. Every derivative of a
multivariate gaussian is written in those two: the score in the mean is
\\w\\, the quadratic form of the log-density is \\r^\top w\\, and the
matrix components are quadratic forms \\w^\top A w\\. Forming them once
per call is what keeps a component from re-solving the same system.

## Usage

``` r
mvg_residuals(y, pc)
```

## Arguments

- y:

  An \\n \times p\\ numeric matrix of observations, already coerced by
  [`as_mv_matrix()`](https://statmodels7.github.io/distributions7/reference/as_mv_matrix.md).

- pc:

  The result of
  [`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md),
  from which `mu` and `sigma_inv` are read.

## Value

A named list with `r` and `w`, each an \\n \times p\\ numeric matrix.
Row \\i\\ of `w` is \\\Sigma^{-1}(y_i - \mu)\\, the right-multiplication
of `r` by `sigma_inv` giving the same rows because \\\Sigma^{-1}\\ is
symmetric.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance and \\y_i\\ the \\i\\-th
row of the response.

## See also

[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
for the argument, and
[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
and
[`distrib_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md)
for the consumers.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
pc <- distributions7:::mvg_pieces(d, theta)
y <- rbind(c(0, 0), c(1, -1), c(2, 0.5))
res <- distributions7:::mvg_residuals(y, pc)

# r is the centered response and w its image under the inverse covariance.
res$r
#>      [,1] [,2]
#> [1,] -0.5  0.3
#> [2,]  0.5 -0.7
#> [3,]  1.5  0.8
all.equal(res$w, res$r %*% solve(pc$sigma), check.attributes = FALSE)
#> [1] TRUE

# The quadratic form of the log-density is the row sum of r * w.
all.equal(rowSums(res$r * res$w),
          mahalanobis(y, pc$mu, pc$sigma))
#> [1] TRUE

# And the score in the mean is w itself, one row per observation.
g <- distrib_gradient(d, y, theta)
all.equal(cbind(g$mu1, g$mu2), res$w, check.attributes = FALSE)
#> [1] TRUE
```
