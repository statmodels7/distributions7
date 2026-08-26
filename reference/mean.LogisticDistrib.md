# Mean of the Logistic Distribution

Closed form: \\E\[Y\] = \mu\\. The density is symmetric about \\\mu\\,
so the location is the mean, the median and the mode at once.

## Arguments

- x:

  A `LogisticDistrib`, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- theta:

  A named list with components `mu` (the location) and `sigma` (the
  scale, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`variance.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.LogisticDistrib.md),
where the scale does enter;
[`kurtosis.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LogisticDistrib.md);
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

## Examples

``` r
d <- logistic_distrib()

# The location is the mean, and the scale does not move it.
mean(d, list(mu = c(-1, 0, 4), sigma = 2))
#> [1] -1  0  4
```
