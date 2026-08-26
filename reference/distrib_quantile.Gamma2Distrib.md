# Gamma Quantile Function in Mean and Variance

Computes the gamma quantile function, the inverse of the regularized
incomplete gamma function in its argument, by calling
[`stats::qgamma()`](https://rdrr.io/r/stats/GammaDist.html) at shape
\\\alpha = \mu^2/\sigma^2\\ and rate \\\lambda = \mu/\sigma^2\\. There
is no elementary closed form;
[`qgamma()`](https://rdrr.io/r/stats/GammaDist.html) inverts the
distribution function numerically. The gamma distribution function is
strictly increasing on \\(0, \infty)\\, so the round trip through
[`distrib_cdf.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma2Distrib.md)
returns `p`.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. Both must be strictly positive.

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
`max(length(p), length(mu), length(sigma2))`.

## See also

[`distrib_cdf.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma2Distrib.md),
which this inverts;
[`distrib_rng.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Gamma2Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
th <- list(mu = 3, sigma2 = 2)

# A central 95 percent interval, visibly asymmetric about the mean of 3.
distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#> [1] 0.9001298 2.7809442 6.3409226

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# The median falls below the mean, the gamma being right skewed.
c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#>   median     mean 
#> 2.780944 3.000000 
```
