# Gaussian Quantile Function in Mean and Precision

Computes the Gaussian quantile function \$\$Q(p; \mu, \tau) = \mu +
\dfrac{\Phi^{-1}(p)}{\sqrt{\tau}}\$\$ by calling
[`stats::qnorm()`](https://rdrr.io/r/stats/Normal.html) at
`sd = 1/sqrt(tau)`. The Gaussian distribution function is strictly
increasing on the whole line, so \\Q\\ is its exact inverse and the
round trip through
[`distrib_cdf.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian3Distrib.md)
returns `p`.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; the endpoints give `-Inf` and `Inf`.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of the length of `p`. A component of length 1 is recycled.
  `tau` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities, which is how a quantile deep in a tail is
  requested without the probability underflowing. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of length
`max(length(p), length(mu), length(tau))`.

## See also

[`distrib_cdf.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian3Distrib.md),
which this inverts;
[`distrib_rng.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Gaussian3Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
th <- list(mu = 1, tau = 0.25)

# The median is mu and the quartiles are symmetric about it.
distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#> [1] -0.3489795  1.0000000  2.3489795

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# Four times the precision halves the distance from the median.
distrib_quantile(d, 0.975, list(mu = 0, tau = 1)) -
  distrib_quantile(d, 0.975, list(mu = 0, tau = 4))
#> [1] 0.979982
```
