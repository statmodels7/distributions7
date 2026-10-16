# Skewness of the Elastic-Net Distribution

Returns 0. The density is symmetric about \\\mu\\ at every parameter
value, so every odd central moment vanishes and the standardized third
one with them. The zero is recycled to the length of the longest
component of `theta`, so the result has one zero per observation when
the components of `theta` vary by observation.

The fourth moment is not zero and is not closed form here:
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
falls to the base class and integrates. Measured at \\\mu = 0\\,
\\\lambda = 2\\ and \\\alpha = 0.5\\ the excess kurtosis is 0.766,
between the Gaussian's 0 and the Laplace's 3.

## Arguments

- x:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- theta:

  A named list with components `mu`, `lambda` and `alpha`. It is aligned
  and validated by name, so a missing or out-of-bounds component signals
  an error.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of the length that the recycled parameters
imply. The family is symmetric about its location at every setting of
the two rates, so the value does not depend on any parameter.

## See also

[`mean.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.EnetDistrib.md)
and
[`variance.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.EnetDistrib.md)
for the other two closed-form moments, and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the fourth.

## Examples

``` r
d <- enet_distrib()
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# The skewness is zero at every parameter value, and a quadrature confirms
# it.
c(closed = skewness(d, th),
  quadrature = integrate(function(u) u^3 * distrib_pdf(d, u, th),
                         -Inf, Inf)$value)
#>     closed quadrature 
#>          0          0 

# The fourth moment is not zero, and sits between the two ends.
c(gaussian = 0, enet = kurtosis(d, th), laplace = 3)
#>  gaussian      enet   laplace 
#> 0.0000000 0.7658603 3.0000000 
```
