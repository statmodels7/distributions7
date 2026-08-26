# Skew Normal Quantile Function in the Centered Parametrization

Computes the quantiles of the skew normal at the direct parameters the
centered triple implies. The skew normal has no closed-form quantile
function, so the value comes from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)'s
root finding on the distribution function, reached through
[`distrib_quantile.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.continuous_distrib.md).

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  when `log.p = TRUE`.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, each a
  numeric vector of length 1 or of the length of `p`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm.
  Defaults to `FALSE`.

- ...:

  Passed to
  [`distrib_quantile.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.continuous_distrib.md),
  including the root finder's tolerance.

## Value

A numeric vector of quantiles, of the length of the recycled inputs.

## See also

[`distrib_cdf.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormal2Distrib.md),
which it inverts,
[`distrib_rng.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormal2Distrib.md)
for draws, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# The round trip through the distribution function.
p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
q <- distrib_quantile(d, p, th)
rbind(quantile = q, back = distrib_cdf(d, q, th))
#>               [,1]       [,2]        [,3]      [,4]     [,5]
#> quantile -1.478544 -0.7072865 -0.09284868 0.6164954 1.792875
#> back      0.050000  0.2500000  0.50000000 0.7500000 0.950000

# A right-skewed density has its median below its mean.
c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#>      median        mean 
#> -0.09284868  0.00000000 
```
