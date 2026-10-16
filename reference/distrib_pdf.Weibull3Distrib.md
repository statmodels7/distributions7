# Weibull Density, Distribution and Generator in the Mean and Shape

The methods of this page describe the Weibull with shape \\\sigma\\ and
scale \\b = m / \Gamma(1 + 1/\sigma)\\, so that \\E\[Y\] = m\\. The
density, the distribution function, the quantile function and the
generator are
[`stats::dweibull()`](https://rdrr.io/r/stats/Weibull.html),
[`stats::pweibull()`](https://rdrr.io/r/stats/Weibull.html),
[`stats::qweibull()`](https://rdrr.io/r/stats/Weibull.html) and
[`stats::rweibull()`](https://rdrr.io/r/stats/Weibull.html) at that
scale. At a small shape \\b\\ underflows to 0; there the four methods
evaluate the same formulas through \\\log b\\ and \\\log(y/b)\\. On the
log scale of the lower tail the distribution function and the quantile
go through the logarithm of the cumulative hazard
([`weibull_log_lower()`](https://statmodels7.github.io/distributions7/reference/weibull_log_lower.md),
[`log_neg_log1mexp()`](https://statmodels7.github.io/distributions7/reference/log_neg_log1mexp.md)),
so that they do not round where the hazard or \\1 - p\\ would.

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
