# Gaussian Log-CDF Hessian in Mean and Variance

Closed form, by the second-order chain rule on the scale
parametrization's derivatives through \\\sigma = \sqrt{\sigma^2}\\. The
map's second partial \\\partial^2\sigma/\partial(\sigma^2)^2 =
-1/(4\sigma^3)\\ contributes to the variance-variance component, which
is the term a first-order chain does not have.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma2` (positive), each a
  numeric vector of length 1 or `n`.

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

\\\mu\\ is the mean, \\\sigma^2 \> 0\\ the variance and \\F\\ the
distribution function.

## See also

[`distrib_grad_cdf.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Gaussian2Distrib.md)
for the first order;
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for the identity;
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

## Examples

``` r
q <- c(-1, 0.5, 2)
th <- list(mu = 0.3, sigma2 = 1.44)

# Against a central difference of the cdf, which shares no arithmetic.
exact <- distrib_hess_cdf(gaussian2_distrib(), q, th, log = FALSE)
fd <- numerical_cdf_deriv(gaussian2_distrib(), q, th, order = 2)
max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#> [1] 6.997102e-08
```
