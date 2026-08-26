# Gaussian Quantile Function

Computes the Gaussian quantile function \$\$Q(p; \mu, \sigma) = \mu +
\sigma\\\Phi^{-1}(p)\$\$ by calling
[`stats::qnorm()`](https://rdrr.io/r/stats/Normal.html). The Gaussian
distribution function is strictly increasing on the whole line, so \\Q\\
is its exact inverse and the round trip through
[`distrib_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian1Distrib.md)
returns `p`.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

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
  logarithms of probabilities, which is how a quantile deep in a tail is
  requested without the probability underflowing. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of length
`max(length(p), length(mu), length(sigma))`.

## See also

[`distrib_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian1Distrib.md),
which this inverts;
[`distrib_rng.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Gaussian1Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- gaussian1_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The median is mu and the quartiles are symmetric about it.
distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#> [1] -0.6117346  0.4000000  1.4117346

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# A quantile 40 standard deviations out, asked for on the log scale
# because the probability itself is zero in double precision.
distrib_quantile(d, pnorm(-40, log.p = TRUE), list(mu = 0, sigma = 1),
                 log.p = TRUE)
#> [1] -40
```
