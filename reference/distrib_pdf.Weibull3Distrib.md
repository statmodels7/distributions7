# Weibull Density, Distribution and Generator in the Mean and Shape

The Weibull with shape \\\sigma\\ and scale \\b = m / \Gamma(1 +
1/\sigma)\\, so that \\E\[Y\] = m\\. The density, the distribution
function, the quantile function and the generator are
[`stats::dweibull()`](https://rdrr.io/r/stats/Weibull.html),
[`stats::pweibull()`](https://rdrr.io/r/stats/Weibull.html),
[`stats::qweibull()`](https://rdrr.io/r/stats/Weibull.html) and
[`stats::rweibull()`](https://rdrr.io/r/stats/Weibull.html) at that
scale.

## Arguments

- distrib:

  A `Weibull3Distrib` object, from
  [`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md).

- y, q:

  A numeric vector of observations or quantiles.

- p:

  A numeric vector of probabilities.

- n:

  The number of draws.

- theta:

  A named list with components `mean` and `sigma`.

- log, log.p, lower.tail:

  As in [`stats::dweibull()`](https://rdrr.io/r/stats/Weibull.html) and
  its siblings.

- ...:

  Unused.

## Value

A numeric vector.

## See also

[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md).

## Examples

``` r
d <- weibull3_distrib()
th <- list(mean = 4, sigma = 1.7)
integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf)$value
#> [1] 4
```
