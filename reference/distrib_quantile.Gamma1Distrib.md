# Gamma Quantile Function in Mean and Dispersion

Computes the gamma quantile function, the inverse of the regularized
incomplete gamma function in its argument, by calling
[`stats::qgamma()`](https://rdrr.io/r/stats/GammaDist.html) at shape \\a
= 1/\phi\\ and rate \\b = 1/(\phi\mu)\\. The gamma distribution function
is strictly increasing on \\(0, \infty)\\, so the round trip through
[`distrib_cdf.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma1Distrib.md)
returns `p`.

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `p`. A component of length 1 is recycled.
  Both must be strictly positive.

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

A numeric vector of quantiles in \\\[0, \infty\]\\, of length
`max(length(p), length(mu), length(phi))`.

## See also

[`distrib_cdf.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma1Distrib.md),
which this inverts;
[`distrib_rng.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Gamma1Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
th <- list(mu = 3, phi = 0.5)

# A central 95 percent interval, visibly asymmetric about the mean of 3.
distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#> [1] 0.3633139 2.5175205 8.3574651

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# The median falls below the mean, the gamma being right skewed.
c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#>  median    mean 
#> 2.51752 3.00000 
```
