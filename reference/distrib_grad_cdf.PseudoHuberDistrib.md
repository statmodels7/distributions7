# Pseudo-Huber Log-CDF Gradient

Closed form in the location and the scale, \\-f(q)\\ and \\-z f(q)\\
with \\z = (q-\mu)/\sigma\\; the shape \\\nu\\ is differenced. The
method is
[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md)
itself, shared with the Student t.

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

A named list of three numeric vectors, `mu`, `sigma` and `nu`, each the
length of `q` recycled against `theta`.

## Details

This family's distribution function is itself a quadrature, so an
evaluation of it is dear and the split is worth more here than
elsewhere: a gradient at 500 quantiles costs 0.08 s against 0.18 s when
all three components are differenced, since the closed pair needs the
density alone.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
shape, \\z = (q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md)
for the shared body;
[`distrib_hess_cdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.PseudoHuberDistrib.md)
for the second order;
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

## Examples

``` r
d <- pseudohuber_distrib()
th <- list(mu = 0.3, sigma = 1.2, nu = 4)
q <- c(-1, 0.5, 2)

# The location component is exact, the density itself.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE

# Differencing the quadrature agrees, and costs more.
fd <- numerical_cdf_deriv(d, q, th, order = 1)
max(abs(fd$mu / distrib_grad_cdf(d, q, th, log = FALSE)$mu - 1))
#> [1] 2.175393e-11
```
