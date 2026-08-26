# Bernoulli Quantile Function

Computes the generalized inverse of the Bernoulli distribution function
by calling [`stats::qbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = 1`. It returns 0 while \\p \le 1 - \mu\\ and 1 afterwards, so
the whole function is one step at \\1 - \mu\\.

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `p`. `mu` must lie in \\(0, 1)\\.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros and ones, of length
`max(length(p), length(mu))`.

## See also

[`distrib_cdf.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BernoulliDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- bernoulli_distrib()
th <- list(mu = 0.3)

# One step, at p = 1 - mu.
distrib_quantile(d, c(0.1, 0.5, 0.69, 0.71, 0.99), th)
#> [1] 0 0 0 1 1

# The median is 0 whenever mu < 1/2 and 1 whenever mu > 1/2.
vapply(c(0.3, 0.7), function(p) distrib_quantile(d, 0.5, list(mu = p)),
       numeric(1))
#> [1] 0 1
```
