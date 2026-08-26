# Logistic Quantile Function

Computes the logistic quantile function \$\$Q(p; \mu, \sigma) = \mu +
\sigma \log\left(\dfrac{p}{1-p}\right)\$\$ by calling
[`stats::qlogis()`](https://rdrr.io/r/stats/Logistic.html). The argument
of the logarithm is the odds, so this function is the logit link applied
to `p` and then carried onto the location and scale.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

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

[`distrib_cdf.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LogisticDistrib.md),
which this inverts;
[`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
which it applies to `p`; and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The median is mu, the family being symmetric.
distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#> [1] -1.247918  0.400000  2.047918

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# It is mu + sigma times the log odds.
all.equal(distrib_quantile(d, p, th), 0.4 + 1.5 * log(p / (1 - p)))
#> [1] TRUE
```
