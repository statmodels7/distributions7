# Chi-Squared Quantile Function

Computes the chi-squared quantile function, the inverse of the
regularized incomplete gamma function in its argument, by calling
[`stats::qchisq()`](https://rdrr.io/r/stats/Chisquare.html) at
`df = mu`. There is no elementary closed form;
[`qchisq()`](https://rdrr.io/r/stats/Chisquare.html) inverts the
distribution function numerically. The distribution function is strictly
increasing on \\(0, \infty)\\, so the round trip through
[`distrib_cdf.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ChisqDistrib.md)
returns `p`.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `p`, recycled if of length 1. It must be strictly
  positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles in \\\[0, \infty\]\\, of length
`max(length(p), length(mu))`.

## See also

[`distrib_cdf.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ChisqDistrib.md),
which this inverts;
[`distrib_rng.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ChisqDistrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
th <- list(mu = 4)

# The critical values a test reads, here at four degrees of freedom.
distrib_quantile(d, c(0.9, 0.95, 0.99), th)
#> [1]  7.779440  9.487729 13.276704

# A central 95 percent interval, asymmetric about the mean of 4.
distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#> [1]  0.4844186  3.3566940 11.1432868

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# The median falls below the mean, the family being right skewed.
c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#>   median     mean 
#> 3.356694 4.000000 
```
