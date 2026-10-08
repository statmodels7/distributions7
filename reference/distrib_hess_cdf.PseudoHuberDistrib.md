# Pseudo-Huber Log-CDF Hessian

Closed form in the location-scale block, and in the mixed components
with a shape \\k\\, \\-f s_k\\ and \\-z f s_k\\ at \\q\\; the components
in the shape parameters alone are integrals of the density's own
derivatives, taken by the compiled rule. The method is
[`compiled_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
itself.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

## Value

A named list of six numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
shape, \\z = (q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`compiled_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
for the shared body;
[`distrib_grad_cdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuberDistrib.md)
for the first order;
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

## Examples

``` r
d <- pseudohuber_distrib()
th <- list(mu = 0.3, sigma = 1.2, nu = 4)

distrib_hess_cdf(d, c(-1, 2), th)$mu_mu
#> [1] -0.1121827 -0.1037853
```
