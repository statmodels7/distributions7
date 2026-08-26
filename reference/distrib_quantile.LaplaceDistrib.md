# Laplace Quantile Function

Computes the Laplace quantile function \$\$Q(p; \mu, \sigma) = \mu -
\sigma\\\mathrm{sign}(p - \tfrac{1}{2})\\\log\left(1 - 2\left\|p -
\tfrac{1}{2}\right\|\right),\$\$ one expression covering both arms. The
median is \\\mu\\, and the two quartiles are \\\mu \pm \sigma \log 2\\.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. The endpoints give `-Inf` and `Inf`; a
  value outside the range gives `NaN`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\ and `p` is replaced by `1 - p` before
  the formula is applied.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are exponentiated
  first. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of length
`max(length(p), length(mu), length(sigma))`.

## See also

[`distrib_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LaplaceDistrib.md),
which this inverts;
[`distrib_rng.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LaplaceDistrib.md),
which draws through it; and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The median is mu and the quartiles are mu -/+ sigma log 2.
distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#> [1] -0.6397208  0.4000000  1.4397208
0.4 + c(-1, 0, 1) * 1.5 * log(2)
#> [1] -0.6397208  0.4000000  1.4397208

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE
```
