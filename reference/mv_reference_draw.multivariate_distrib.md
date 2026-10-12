# An Inflated Gaussian Proposal

The default proposal: a gaussian with the family's mean and twice its
covariance, drawn as \\\mu + Lz\\ with \\LL^\top = 2\Sigma\\. For a
gaussian family a proposal with the family's own covariance would make
every importance ratio one, and the inflation keeps the ratio
informative.

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

The covariance must be finite and non-singular, which a family supported
on all of \\\mathbb{R}^p\\ with finite second moments has. A family
whose support is a lower-dimensional set, or whose covariance is
infinite, registers its own method (see
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)).

The log-density is computed from the Cholesky factor \\L\\: the whitened
residual \\L^{-1}r\\ solves a lower triangular system, which
[`base::forwardsolve()`](https://rdrr.io/r/base/backsolve.html)
evaluates.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance, \\L\\ a lower Cholesky
factor of \\2\Sigma\\ and \\p\\ the dimension.

## See also

[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
for the generic and the proposals of the families that register their
own, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
for the consumer.

## Examples

``` r
d <- mvgaussian1_distrib(2)
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

# logd is the log-density of the inflated gaussian at those draws.
S2 <- 2 * variance(d, theta)
mu <- c(1, -1)
ref <- -0.5 * (2 * log(2 * pi) + log(det(S2)) + mahalanobis(r$y, mu, S2))
all.equal(r$logd, as.numeric(ref))
#> [1] TRUE

# The sample is twice as spread as the family.
set.seed(2)
round(rbind(proposal = diag(var(mv_reference_draw(d, theta, 20000)$y)),
            family = diag(variance(d, theta))), 3)
#>             v1    v2
#> proposal 2.024 2.526
#> family   1.000 1.250
```
