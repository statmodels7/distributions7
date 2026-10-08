# Student t Density, Distribution and Generator in the Standard Deviation

The Student t with \\\nu\\ degrees of freedom, location \\\mu\\ and
scale \\s_0 = \sigma\sqrt{1 - 2/\nu}\\, so that the standard deviation
is \\\sigma\\. The density, the distribution function, the quantile
function and the generator are
[`stats::dt()`](https://rdrr.io/r/stats/TDist.html),
[`stats::pt()`](https://rdrr.io/r/stats/TDist.html),
[`stats::qt()`](https://rdrr.io/r/stats/TDist.html) and
[`stats::rt()`](https://rdrr.io/r/stats/TDist.html) at that scale.

## Arguments

- distrib:

  A `StudentT2Distrib` object, from
  [`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md).

- y, q:

  A numeric vector of observations or quantiles.

- p:

  A numeric vector of probabilities.

- n:

  The number of draws.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- log, log.p, lower.tail:

  As in [`stats::dt()`](https://rdrr.io/r/stats/TDist.html) and its
  siblings.

- ...:

  Unused.

## Value

A numeric vector.

## See also

[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md).

## Examples

``` r
d <- student_t2_distrib()
th <- list(mu = 1, sigma = 2, nu = 6)
integrate(function(t) (t - 1)^2 * distrib_pdf(d, t, th), -Inf, Inf)$value
#> [1] 4
```
