# Lognormal Random Number Generator

Draws `n` independent lognormal variates by calling
[`stats::rlnorm()`](https://rdrr.io/r/stats/Lognormal.html) at
`meanlog = mu` and `sdlog = sqrt(sigma2)`, so the draws come from R's
own generator and depend on `.Random.seed` in the usual way. The
ratio-of-uniforms fallback the base class supplies is bypassed.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled, so
  a vector of length `n` draws one variate per parameter setting.
  `sigma2` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` strictly positive draws.

## See also

[`distrib_quantile.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Lognormal1Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()

# Same generator as stats::rlnorm at this parametrization.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 0.5, sigma2 = 0.36))
set.seed(2)
identical(a, rlnorm(3, 0.5, sqrt(0.36)))
#> [1] TRUE

# The parameters are moments of the LOGARITHM, so that is where the sample
# recovers them.
set.seed(6)
z <- distrib_rng(d, 2e4, list(mu = 0.5, sigma2 = 0.36))
c(mu = mean(log(z)), sigma2 = var(log(z)))
#>        mu    sigma2 
#> 0.4990354 0.3598171 

# On the original scale the sample mean is exp(mu + sigma2/2).
c(sample = mean(z), theory = exp(0.5 + 0.36 / 2))
#>   sample   theory 
#> 1.972276 1.973878 
```
