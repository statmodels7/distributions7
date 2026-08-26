# Gaussian Log-CDF Hessian in Mean and Precision

Closed form, by the second-order chain rule through \\\sigma =
\tau^{-1/2}\\. The map's second partial
\\\partial^2\sigma/\partial\tau^2 = \tfrac34\tau^{-5/2}\\ contributes to
the precision-precision component.

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

A named list of three numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the mean, \\\tau \> 0\\ the precision and \\F\\ the
distribution function.

## See also

[`distrib_grad_cdf.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Gaussian3Distrib.md)
for the first order;
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for the identity;
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

## Examples

``` r
q <- c(-1, 0.5, 2)
th <- list(mu = 0.3, tau = 1 / 1.44)

# Against a central difference of the cdf, which shares no arithmetic.
exact <- distrib_hess_cdf(gaussian3_distrib(), q, th, log = FALSE)
fd <- numerical_cdf_deriv(gaussian3_distrib(), q, th, order = 2)
max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#> [1] 7.91399e-08
```
