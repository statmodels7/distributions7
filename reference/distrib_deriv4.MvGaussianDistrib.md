# Multivariate Gaussian Fourth Derivatives

Computes every fourth derivative of the log-density in the parameters,
in closed form, by the same algebra as
[`distrib_deriv3.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvGaussianDistrib.md)
one order up: the tuple is split into mean indices and matrix indices,
three or more mean indices give exactly zero, and the rest are read off
the precision's derivative array \\P_t\\ and the fourth derivative of
the log-determinant. The arrays come from
[`parameters7::param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.html)
under a precision parametrization and from the expansion of the
derivative of an inverse under a covariance one.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. With
  `expected = TRUE` only its row count is used.

- theta:

  A named list of parameters, each component a single number.

- expected:

  Logical of length 1. When `TRUE` the expectation of each component is
  returned, by sampling. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. Every link of this family is the identity, so the two
  coincide.

- approx:

  Ignored: sampling is the only multivariate route to an expectation
  here. Present so that the signature matches the generic's.

- nsim:

  The number of draws used when `expected = TRUE`. Defaults to `10000`.
  Ignored otherwise.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 4)`. With `expected = TRUE` every vector is
constant.

## Details

With `expected = TRUE` the expectation is taken by sampling `nsim` draws
from the family and averaging the observed components, exactly as at
order three. `approx` is not read, and the result carries Monte Carlo
error of order `nsim^(-1/2)`.

## Notation

\\\mu\\ is the mean, \\M\\ the matrix the parametrization carries,
\\\eta\\ its free vector, \\r = y - \mu\\ the centered response and
\\P_t\\ the precision's derivative array over the multiset \\t\\.

## See also

[`distrib_deriv3.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvGaussianDistrib.md)
for the order below,
[`mvg_higher()`](https://statmodels7.github.io/distributions7/reference/mvg_higher.md)
for the shared engine, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)

d4 <- distrib_deriv4(d, y, theta)
length(d4)
#> [1] 70

# Four mean indices vanish, as three already did.
d4[["mu1_mu1_mu2_mu2"]]
#> [1] 0 0 0 0

# Against one stencil on the analytic third order.
h <- 1e-4
tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
c(exact = sum(d4[["mu1_mu2_sigma_L2.1_sigma_L2.1"]]),
  stencil = (sum(distrib_deriv3(d, y, tp)[["mu1_mu2_sigma_L2.1"]]) -
             sum(distrib_deriv3(d, y, tm)[["mu1_mu2_sigma_L2.1"]])) / (2 * h))
#>        exact      stencil 
#> 1.776357e-15 4.440892e-12 
```
