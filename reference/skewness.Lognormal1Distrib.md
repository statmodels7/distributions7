# Skewness of the Lognormal Distribution

Closed form: \\\gamma_1 = (e^{\sigma^2} + 2)\sqrt{e^{\sigma^2} - 1}\\,
free of \\\mu\\. It is positive at every parameter value and grows
without bound with the spread of the logarithm: at \\\sigma^2 = 1\\ it
is already 6.18, where a gamma of the same mean would need a shape near
0.1 to match it.

## Arguments

- x:

  A `Lognormal1Distrib`, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- theta:

  A named list with components `mu` (any real value) and `sigma2`
  (positive), each a numeric vector of length 1 or `n`. Only `sigma2`
  enters the value, and the exponentials overflow to `Inf` for large
  `sigma2`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$mu), length(theta$sigma2))`, positive throughout.

## Notation

\\\sigma^2 \> 0\\ is the variance of \\\log Y\\. The location of the
logarithm does not enter a standardized moment.

## See also

[`kurtosis.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Lognormal1Distrib.md),
which grows faster still;
[`variance.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Lognormal1Distrib.md);
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

## Examples

``` r
d <- lognormal1_distrib()

# The published form, written out.
all.equal(skewness(d, list(mu = 0, sigma2 = 1)),
          (exp(1) + 2) * sqrt(exp(1) - 1))
#> [1] TRUE

# The location of the logarithm does not enter it.
skewness(d, list(mu = c(-2, 0, 5), sigma2 = 0.25))
#> [1] 1.75019 1.75019 1.75019
```
