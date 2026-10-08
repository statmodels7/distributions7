# Generalized Gamma Density, Distribution and Generator in the Mean

The generalized gamma with shape \\d\\, power \\p\\ and scale \\a =
m\\\Gamma(d/p)/\Gamma((d+1)/p)\\, so that \\E\[Y\] = m\\. The density,
the distribution function, the quantile function and the generator are
those of
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
at that scale: \\(Y/a)^p\\ is a gamma variable with shape \\d/p\\ and
unit rate.

## Arguments

- distrib:

  A `GenGamma2Distrib` object, from
  [`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md).

- y, q:

  A numeric vector of observations or quantiles.

- p:

  A numeric vector of probabilities.

- n:

  The number of draws.

- theta:

  A named list with components `mean`, `d` and `p`.

- log, log.p, lower.tail:

  As in [`stats::pgamma()`](https://rdrr.io/r/stats/GammaDist.html) and
  its siblings.

- ...:

  Unused.

## Value

A numeric vector.

## See also

[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md).

## Examples

``` r
d <- gengamma2_distrib()
th <- list(mean = 5, d = 3, p = 1.5)
integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf)$value
#> [1] 5
```
