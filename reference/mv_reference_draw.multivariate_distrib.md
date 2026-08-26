# An Inflated Gaussian Proposal

The default proposal: a gaussian with the family's mean and TWICE its
covariance, drawn as \\\mu + Lz\\ with \\LL^\top = 2\Sigma\\. The
inflation makes the importance ratio informative; at the family's own
covariance every ratio would be one.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with a mean and a finite covariance.

- theta:

  A named list of parameters, already aligned by the generic.

- n:

  The number of draws, a single positive whole number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `y`, an \\n \times p\\ numeric matrix, and `logd`, a
numeric vector of length \\n\\.

## Details

It requires the covariance to be non-singular, which a support of full
dimension in \\\mathbb{R}^p\\ gives. On a family whose support is a
lower-dimensional set the covariance is singular and
[`chol()`](https://rdrr.io/r/base/chol.html) may accept it anyway, so
the failure is quiet: see
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
for the measured consequence. Such a family registers its own method.

The log-density is formed from the Cholesky factor directly. The
whitened residual is \\L^{-1}r\\, a LOWER triangular system, so
[`base::forwardsolve()`](https://rdrr.io/r/base/backsolve.html) is what
solves it; [`backsolve()`](https://rdrr.io/r/base/backsolve.html) on the
transpose solves \\L^\top x = b\\, a different vector of the same shape.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance, \\L\\ a lower Cholesky
factor of \\2\Sigma\\ and \\p\\ the dimension.

## See also

[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
for the generic and the override a simplex-valued family needs, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
for the consumer.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)

set.seed(1)
r <- mv_reference_draw(d, theta, 4)
r$y
#>              v1         v2
#> [1,]  0.1140605 -0.9769754
#> [2,]  1.2597109 -2.0304621
#> [3,] -0.1817573 -0.9015499
#> [4,]  3.2560677  1.1721827

# logd really is the log-density of the inflated gaussian at those draws.
S2 <- 2 * variance(d, theta)
mu <- c(1, -1)
ref <- -0.5 * (2 * log(2 * pi) + log(det(S2)) + mahalanobis(r$y, mu, S2))
all.equal(r$logd, as.numeric(ref))
#> [1] TRUE

# And the sample is twice as spread as the family, which is the point.
set.seed(2)
round(rbind(proposal = diag(var(mv_reference_draw(d, theta, 20000)$y)),
            family = diag(variance(d, theta))), 3)
#>             v1    v2
#> proposal 2.024 2.526
#> family   1.000 1.250
```
