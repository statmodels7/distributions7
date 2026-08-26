# Multivariate Gaussian Density

Computes the multivariate gaussian log-density \$\$\ell =
-\frac{p}{2}\log 2\pi - \frac{1}{2}\log\lvert\Sigma\rvert -
\frac{1}{2}(y-\mu)^\top \Sigma^{-1} (y-\mu),\$\$ row by row, and
exponentiates it unless `log = TRUE`. The log-determinant comes from the
matrix parametrization's own `param_logdet()` and the quadratic form
from its `param_solve()`, so neither a determinant nor an explicit
inverse is formed by this method. Both are computed once per call, not
once per observation, since the matrix does not vary with the response.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations, one row each. A
  plain numeric vector of length \\p\\ is read as a single observation.
  Every point of \\\mathbb{R}^p\\ is in the support, so no row is
  rejected. A matrix with zero rows returns `numeric(0)`.

- theta:

  A named list of parameters, each component a single number: `mu1`,
  ..., `mup` and the matrix parametrization's prefixed free values. A
  parameter may not vary by observation here, and a component longer
  than one is an error.

- log:

  Logical of length 1. When `TRUE` the log-density is returned, which
  stays finite for a point the density itself underflows at. Defaults to
  `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length \\n\\, one density per row of `y`.

## See also

[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
for the derivatives of this log-density,
[`distrib_rng.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvGaussianDistrib.md)
to draw from it,
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md)
for the covariance it uses, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
y <- rbind(c(0, 0), c(1, -1), c(0.5, -0.3))

distrib_pdf(d, y, theta)
#> [1] 0.13361803 0.08899957 0.17589341

# Against the formula written out with an explicit inverse.
S <- mv_sigma(d, theta)
mu <- c(0.5, -0.3)
ref <- -0.5 * (2 * log(2 * pi) + log(det(S)) + mahalanobis(y, mu, S))
all.equal(distrib_pdf(d, y, theta, log = TRUE), as.numeric(ref))
#> [1] TRUE

# A vector of length p is one observation, and the mode is the mean.
distrib_pdf(d, c(0.5, -0.3), theta, log = TRUE)
#> [1] -1.737877

# Far out in the tail the density underflows and its logarithm does not.
distrib_pdf(d, c(60, -60), theta)
#> [1] 0
distrib_pdf(d, c(60, -60), theta, log = TRUE)
#> [1] -6373.378
```
