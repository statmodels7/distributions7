# Pseudo-Huber Log-CDF Gradient

Closed form in the location and the scale, \\-f(q)\\ and \\-z f(q)\\
with \\z = (q-\mu)/\sigma\\; the component in the shape \\\nu\\ is the
integral of the density's own derivative, taken by the compiled rule.
The method is
[`compiled_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
itself, shared with the other families whose distribution function has
no closed derivative in a shape parameter.

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

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
shape, \\z = (q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`compiled_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
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

# The shape component, on the log scale of the lower tail.
distrib_grad_cdf(d, q, th)$nu
#> [1]  0.033358345 -0.003435126 -0.011536871
```
