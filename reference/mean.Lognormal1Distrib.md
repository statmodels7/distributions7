# Mean of the Lognormal Distribution

Closed form: \\E\[Y\] = \exp(\mu + \sigma^2/2)\\. The parameters
describe \\\log Y\\ and not \\Y\\, so \\\mu\\ is the mean of the
logarithm and \\\exp(\mu)\\ is the median of \\Y\\; the mean is larger
than the median by the factor \\\exp(\sigma^2/2)\\, which grows with the
spread.

## Arguments

- x:

  A `Lognormal1Distrib`, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- theta:

  A named list with components `mu` (the mean of the logarithm, any real
  value) and `sigma2` (its variance, positive), each a numeric vector of
  length 1 or `n`. The mean overflows to `Inf` for large `sigma2`, the
  exponential of half of it being taken.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma2))`.

## Notation

\\\mu\\ is the mean of \\\log Y\\ and \\\sigma^2 \> 0\\ its variance.

## See also

[`variance.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Lognormal1Distrib.md),
[`skewness.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Lognormal1Distrib.md);
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the median, which is \\\exp(\mu)\\;
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

## Examples

``` r
d <- lognormal1_distrib()

# exp(mu + sigma2 / 2), which exceeds the median exp(mu).
all.equal(mean(d, list(mu = 0, sigma2 = 1)), exp(0.5))
#> [1] TRUE
distrib_quantile(d, 0.5, list(mu = 0, sigma2 = 1))
#> [1] 1
```
