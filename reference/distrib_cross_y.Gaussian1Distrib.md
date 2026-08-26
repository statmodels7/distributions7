# Gaussian Mixed Derivatives

Closed form. With \\r = y - \mu\\, \$\$\frac{\partial^2 \ell}{\partial
y\\ \partial \mu} = \frac{1}{\sigma^2}, \qquad \frac{\partial^2
\ell}{\partial y\\ \partial \sigma} = \frac{2r}{\sigma^3}.\$\$

The first is constant along `y`, a quadratic log-density in a location
parameter giving nothing else; the second changes sign with the
residual, being zero at the mean and growing linearly away from it.

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, aligned by the generic.
  Each may be length 1 or `length(y)`.

- scale:

  Handled by the generic after dispatch; this method always returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma`, each a numeric vector of
length `length(y)`.

## See also

[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic;
[`distrib_cross_y.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.StudentT1Distrib.md),
whose mu component is not constant;
[`distrib_grad_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gaussian1Distrib.md)
for the quantity differentiated.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
distrib_cross_y(d, y, list(mu = 0.3, sigma = 1.4))
#> $mu
#> [1] 0.5102041 0.5102041 0.5102041
#> 
#> $sigma
#> [1] -0.9475219 -0.2186589  1.2390671
#> 

# The sigma component vanishes exactly at the mean.
distrib_cross_y(d, 0.3, list(mu = 0.3, sigma = 1.4))$sigma
#> [1] 0
```
