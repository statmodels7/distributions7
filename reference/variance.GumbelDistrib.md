# Variance of the Gumbel Distribution

Closed form: \\\operatorname{Var}(Y) = \pi^2\sigma^2/6\\. The scale is
not the standard deviation: it is smaller by \\\pi/\sqrt6 \approx
1.2825\\. The location does not enter, the family being location-scale.

## Arguments

- x:

  A `GumbelDistrib`, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- theta:

  A named list with components `mu` (the location) and `sigma` (the
  scale, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`. The location does not
enter the value, so a setting that varies `mu` alone repeats one number.

## See also

[`mean.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.GumbelDistrib.md),
[`kurtosis.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GumbelDistrib.md),
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

## Examples

``` r
d <- gumbel_distrib()

# pi^2 sigma^2 / 6.
all.equal(variance(d, list(mu = 0, sigma = 2)), pi^2 * 4 / 6)
#> [1] TRUE

# The standard deviation is pi / sqrt(6) times the scale.
std_dev(d, list(mu = 0, sigma = 1))
#> [1] 1.28255
```
