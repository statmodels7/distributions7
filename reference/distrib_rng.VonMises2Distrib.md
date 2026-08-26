# von Mises Random Generation in the Resultant Length

Draws `n` independent angles by inverting the map once and calling
[`distrib_rng.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.VonMises1Distrib.md)
at the implied concentration \\\kappa = A^{-1}(\rho)\\. The generator is
the rejection algorithm of Best and Fisher (1979), which involves no
Bessel function; the only Bessel work here is the single inversion of
\\A\\.

## Arguments

- distrib:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1. `mu` must lie in \\(-\pi, \pi)\\ and `rho` in \\(0, 1)\\.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` angles in \\\[-\pi, \pi)\\.

## See also

[`distrib_rng.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.VonMises1Distrib.md),
which this calls;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample; and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
th <- list(mu = 0.5, rho = 0.7)
set.seed(1)
z <- distrib_rng(d2, 3e5, th)

# This parametrization is the one a sample reads back directly: the mean
# resultant length of the draws is rho, and the circular mean is mu.
c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2), rho = 0.7)
#> resultant       rho 
#> 0.7003239 0.7000000 
c(circular_mean = atan2(mean(sin(z)), mean(cos(z))), mu = 0.5)
#> circular_mean            mu 
#>     0.4990351     0.5000000 
```
