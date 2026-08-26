# Mean of the Gumbel Distribution

Closed form: \\E\[Y\] = \mu + \gamma\sigma\\, with \\\gamma \approx
0.5772157\\ the Euler-Mascheroni constant. The location is therefore not
the mean: the distribution of a maximum is shifted to the right of its
location by a fixed multiple of the scale. The constant is obtained as
`-digamma(1)`, which is exact to the last bit.

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

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale and \\\gamma\\ the
Euler-Mascheroni constant.

## See also

[`variance.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.GumbelDistrib.md);
[`skewness.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GumbelDistrib.md)
and
[`kurtosis.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GumbelDistrib.md),
which are constants;
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

## Examples

``` r
d <- gumbel_distrib()

# The location plus the Euler-Mascheroni constant times the scale.
all.equal(mean(d, list(mu = 0, sigma = 1)), -digamma(1))
#> [1] TRUE

# The gap between location and mean grows with the scale.
mean(d, list(mu = 0, sigma = c(1, 2, 5)))
#> [1] 0.5772157 1.1544313 2.8860783
```
