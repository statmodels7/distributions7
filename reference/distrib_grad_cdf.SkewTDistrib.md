# Skew t Log-CDF Gradient

Closed form in the location and the scale, \\-f(q)\\ and \\-z f(q)\\
with \\z = (q-\mu)/\sigma\\; the component in the shape and in the
degrees of freedom are integrals, each is the integral of the density's
own derivative, taken by the compiled rule. The method is
[`compiled_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
itself, shared with the other families whose distribution function has
no closed derivative in a shape parameter.

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

A named list of four numeric vectors, `mu`, `sigma`, `alpha` and `nu`,
each the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\alpha\\ the
shape, \\\nu \> 0\\ the degrees of freedom, \\z = (q-\mu)/\sigma\\ and
\\f\\ the density.

## See also

[`compiled_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
for the shared body;
[`distrib_hess_cdf.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.SkewTDistrib.md)
for the second order;
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

## Examples

``` r
d <- skewt_distrib()
th <- list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6)
q <- c(-1, 0.5, 2)

# The location component is exact, the density itself.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE
```
