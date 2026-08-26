# Gaussian Log-CDF Gradient in Mean and Precision

Closed form, by the chain rule on the scale parametrization's
derivatives through \\\sigma = \tau^{-1/2}\\. The mean component is
unchanged, \\-f(q)\\; the precision component is the scale one times
\\-\tfrac12\tau^{-3/2}\\, so it carries the opposite sign, a larger
precision being a tighter distribution.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `tau` (positive), each a numeric
  vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `tau`, each the length of
`q` recycled against `theta`.

## Notation

\\\mu\\ is the mean, \\\tau \> 0\\ the precision, \\f\\ the density and
\\F\\ the distribution function.

## See also

[`distrib_hess_cdf.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.Gaussian3Distrib.md)
for the second order;
[`distrib_grad_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Gaussian1Distrib.md),
the parent;
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

## Examples

``` r
q <- c(-1, 0.5, 2)
g3 <- distrib_grad_cdf(gaussian3_distrib(), q,
                       list(mu = 0.3, tau = 1 / 1.44), log = FALSE)
g1 <- distrib_grad_cdf(gaussian1_distrib(), q,
                       list(mu = 0.3, sigma = 1.2), log = FALSE)

# The precision component is the scale one times -0.5 tau^(-3/2).
all.equal(g3$tau, g1$sigma * (-0.5) * (1 / 1.44)^(-1.5))
#> [1] TRUE
```
