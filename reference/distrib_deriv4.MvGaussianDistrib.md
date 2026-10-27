# Multivariate Gaussian Fourth Derivatives

Computes every fourth derivative of the log-density in the parameters,
in closed form, by the algebra of
[`distrib_deriv3.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MvGaussianDistrib.md)
one order up: the tuple is split into mean indices and matrix indices,
three or more mean indices give exactly zero, and the rest are read off
the derivative of the precision \\P_t\\ from
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md)
and the fourth derivative of the log-determinant.

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

  A named list of parameters in which each component is a single number.

- expected:

  Logical of length 1. When `TRUE` the expectation of each component is
  returned. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`. The method returns the
  parameter scale and the generic applies the link scale to its result.
  Every link of this family is the identity, so the two coincide.

- method:

  `NULL`: the generic applies an object before dispatching on the
  family, so this method receives `NULL` and returns the closed form.
  See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 4)`. With `expected = TRUE` every vector is
constant.

## Details

With `expected = TRUE` the expectation is exact, for the reason given at
order three.

## Notation

\\\mu\\ is the mean, \\M\\ the matrix that the parametrization carries,
\\\eta\\ its free vector, \\r = y - \mu\\ the centered response and
\\P_t\\ the derivative of the precision over the multiset \\t\\.

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

# The result is compared with one stencil on the analytic third order.
h <- 1e-4
tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
c(exact = sum(d4[["mu1_mu2_sigma_L2.1_sigma_L2.1"]]),
  stencil = (sum(distrib_deriv3(d, y, tp)[["mu1_mu2_sigma_L2.1"]]) -
             sum(distrib_deriv3(d, y, tm)[["mu1_mu2_sigma_L2.1"]])) / (2 * h))
#>   exact stencil 
#>       0       0 
```
