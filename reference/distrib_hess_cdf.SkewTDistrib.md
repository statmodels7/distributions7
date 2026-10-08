# Skew t Log-CDF Hessian

Closed form in the location-scale block, and in the mixed components
with a shape \\k\\, \\-f s_k\\ and \\-z f s_k\\ at \\q\\; the components
in the shape parameters alone are integrals of the density's own
derivatives, taken by the compiled rule. The method is
[`compiled_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
itself.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive), `alpha` (any
  sign) and `nu` (positive), each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

## Value

A named list of ten numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\alpha\\ the
shape, \\\nu \> 0\\ the degrees of freedom, \\z = (q-\mu)/\sigma\\ and
\\f\\ the density.

## See also

[`compiled_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
for the shared body;
[`distrib_grad_cdf.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.SkewTDistrib.md)
for the first order;
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

## Examples

``` r
d <- skewt_distrib()
th <- list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6)

# Ten components, keyed as hess_names().
names(distrib_hess_cdf(d, c(-1, 2), th))
#>  [1] "mu_mu"       "sigma_sigma" "alpha_alpha" "nu_nu"       "mu_sigma"   
#>  [6] "mu_alpha"    "mu_nu"       "sigma_alpha" "sigma_nu"    "alpha_nu"   
```
