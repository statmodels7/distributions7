# Beta Quantile Function in the Shapes

Computes the beta quantile function, the inverse of the regularized
incomplete beta function in its argument, by calling
[`stats::qbeta()`](https://rdrr.io/r/stats/Beta.html) at
`shape1 = alpha` and `shape2 = beta`. There is no elementary closed
form; [`qbeta()`](https://rdrr.io/r/stats/Beta.html) inverts the
distribution function numerically. The beta distribution function is
strictly increasing on \\(0, 1)\\, so the round trip through
[`distrib_cdf.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Beta2Distrib.md)
returns `p`.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; the endpoints give 0 and 1.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. Both must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles in \\\[0, 1\]\\, of length
`max(length(p), length(alpha), length(beta))`.

## See also

[`distrib_cdf.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Beta2Distrib.md),
which this inverts;
[`distrib_rng.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Beta2Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
th <- list(alpha = 2, beta = 5)

# A central 95 percent interval, asymmetric about the mean of 2/7.
distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#> [1] 0.04327187 0.26444998 0.64123458

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# Both shapes 1 is the uniform, so Q(p) = p.
distrib_quantile(d, p, list(alpha = 1, beta = 1))
#> [1] 0.025 0.500 0.975
```
