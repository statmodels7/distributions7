# Cauchy Quantile Function

Computes the Cauchy quantile function \$\$Q(p; \mu, \sigma) = \mu +
\sigma \tan\left(\pi\left(p - \dfrac{1}{2}\right)\right)\$\$ by calling
[`stats::qcauchy()`](https://rdrr.io/r/stats/Cauchy.html). The median is
\\\mu\\ and the quartiles are \\\mu \pm \sigma\\, so the two parameters
are read off this function directly: they are the median and the
half-interquartile range, the location and spread this family has in
place of a mean and a standard deviation.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; the endpoints give `-Inf` and `Inf`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of length
`max(length(p), length(mu), length(sigma))`.

## See also

[`distrib_cdf.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.CauchyDistrib.md),
which this inverts;
[`skewness.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.CauchyDistrib.md)
for why the median and the half-interquartile range are what this family
offers;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The median is mu and the quartiles are mu -/+ sigma exactly.
distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#> [1] -1.1  0.4  1.9

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.01, 0.5, 0.99)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# The tails are heavy, so the extreme quantiles are far out: the 99.9th
# percentile sits at over 300 scale units.
(distrib_quantile(d, 0.999, th) - 0.4) / 1.5
#> [1] 318.3088
```
