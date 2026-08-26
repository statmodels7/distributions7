# Laplace Quantile Function, Rate Parametrization

Computes the Laplace quantile function in the rate parametrization,
\$\$Q(p; \mu, \lambda) = \mu - \dfrac{\mathrm{sign}(p -
\tfrac{1}{2})}{\lambda}\\\log\left(1 - 2\left\|p -
\tfrac{1}{2}\right\|\right),\$\$ one expression covering both arms. The
median is \\\mu\\ and the quartiles are \\\mu \pm \log(2)/\lambda\\.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. The endpoints give `-Inf` and `Inf`.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. `lambda` must be strictly positive.

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
`max(length(p), length(mu), length(lambda))`.

## See also

[`distrib_cdf.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Laplace2Distrib.md),
which this inverts;
[`distrib_rng.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Laplace2Distrib.md),
which draws through it; and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
th <- list(mu = 0.4, lambda = 2)

# The median is mu and the quartiles are mu -/+ log(2)/lambda.
distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#> [1] 0.05342641 0.40000000 0.74657359
0.4 + c(-1, 0, 1) * log(2) / 2
#> [1] 0.05342641 0.40000000 0.74657359

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE
```
