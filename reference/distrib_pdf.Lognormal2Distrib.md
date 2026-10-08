# Lognormal Density, Distribution and Generator in the Mean and Variance of Y

The lognormal with \\\log Y \sim N(\mu_l, S)\\, where \\S = \log(1 +
v/m^2)\\ and \\\mu_l = \log m - S/2\\, so that \\E\[Y\] = m\\ and
\\\operatorname{Var}(Y) = v\\. The density, the distribution function,
the quantile function and the generator are
[`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html),
[`stats::plnorm()`](https://rdrr.io/r/stats/Lognormal.html),
[`stats::qlnorm()`](https://rdrr.io/r/stats/Lognormal.html) and
[`stats::rlnorm()`](https://rdrr.io/r/stats/Lognormal.html) at `meanlog`
\\\mu_l\\ and `sdlog` \\\sqrt{S}\\, with \\S\\ formed by
[`log1p()`](https://rdrr.io/r/base/Log.html).

## Arguments

- distrib:

  A `Lognormal2Distrib` object, from
  [`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md).

- y, q:

  A numeric vector of observations or quantiles.

- p:

  A numeric vector of probabilities.

- n:

  The number of draws.

- theta:

  A named list with components `mean` and `var`.

- log, log.p, lower.tail:

  As in [`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html) and
  its siblings.

- ...:

  Unused.

## Value

A numeric vector.

## See also

[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md).

## Examples

``` r
d <- lognormal2_distrib()
th <- list(mean = 3, var = 2)
integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf)$value
#> [1] 3
distrib_quantile(d, distrib_cdf(d, c(1, 3), th), th)
#> [1] 1 3
```
