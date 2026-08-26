# Generalized Pareto Density

Computes the generalized Pareto density \$\$f(y; \sigma, \xi) =
\dfrac{1}{\sigma} \left(1 + \dfrac{\xi y}{\sigma}\right)^{-1/\xi - 1},
\qquad y \ge 0,\\ 1 + \xi y/\sigma \> 0,\$\$ and 0 outside that region.
At \\\xi = 0\\ the limit is the exponential density
\\e^{-y/\sigma}/\sigma\\.

The compiled kernel does not branch on \\\xi = 0\\. It writes the
log-survival as \\-(y/\sigma)\Lambda(u)\\ with \\u = \xi y/\sigma\\ and
\\\Lambda(u) = \log(1+u)/u\\, which is analytic with \\\Lambda(0) = 1\\,
so every division by the shape disappears and the exponential limit is
an ordinary point of the formula. Measured, the density agrees with
`dexp` to \\1.4\times10^{-17}\\ at \\\xi = 0\\ exactly.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector of observations. A negative value, or one beyond
  [`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md),
  gives a density of 0.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`. `sigma` must be strictly
  positive; `xi` may take either sign.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A numeric vector of densities, of the length of the recycled inputs.

## Notation

\\\sigma \> 0\\ is the scale and \\\xi\\ the shape; \\\sigma\\ is not
the mean, which is \\\sigma/(1-\xi)\\ and exists only for \\\xi \< 1\\.

## See also

[`distrib_cdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GPDDistrib.md)
for the distribution function,
[`distrib_gradient.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GPDDistrib.md)
for the score,
[`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md)
for where the support ends, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
y <- c(0.2, 1, 4)
th <- list(sigma = 1.5, xi = 0.3)

# The formula written out.
all.equal(distrib_pdf(d, y, th), (1 + 0.3 * y / 1.5)^(-1 / 0.3 - 1) / 1.5)
#> [1] TRUE

# It integrates to one.
integrate(function(v) distrib_pdf(d, v, th), 0, Inf)$value
#> [1] 1

# Shape zero is the exponential, and the limit is reached by a series
# rather than by a branch: the error is linear in xi.
vapply(c(0, 1e-10, 1e-6, 1e-3), function(x)
  max(abs(distrib_pdf(d, y, list(sigma = 1.5, xi = x)) -
          dexp(y, rate = 1 / 1.5))), 0)
#> [1] 1.387779e-17 1.521239e-11 1.521235e-07 1.520476e-04

# A negative shape bounds the support at -sigma/xi.
distrib_pdf(d, c(4, 5, 6), list(sigma = 2, xi = -0.4))
#> [1] 0.04472136 0.00000000 0.00000000
```
