# Weibull First Derivative in the Response

Computes \\\partial \ell / \partial y\\, the derivative of the Weibull
log-density with respect to the response, in closed form. With \\u =
(y/\mu)^{\sigma}\\, \$\$\dfrac{\partial \ell}{\partial y} =
\dfrac{\sigma - 1 - \sigma u}{y}.\$\$ At \\\sigma \> 1\\ it vanishes at
the mode \\y = \mu\\(\sigma-1)/\sigma\\^{1/\sigma}\\; at \\\sigma \le
1\\ the density is decreasing on the whole support and the derivative is
negative everywhere. This quantity is what a quantile residual's
delta-method standard error and a change of variable in the response
both need.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- y:

  A numeric vector of positive observations. At or below zero the result
  is not defined and propagates as `NaN` or `Inf`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(sigma))`,
one value per observation.

## See also

[`distrib_hess_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Weibull1Distrib.md)
for the second derivative in the response,
[`distrib_cross_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Weibull1Distrib.md)
for the mixed derivative in the response and the parameters,
[`distrib_gradient.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull1Distrib.md)
for the derivatives in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)

# The closed form, written out.
u <- (y / 2)^1.5
all.equal(distrib_grad_y(d, y, th), (1.5 - 1 - 1.5 * u) / y)
#> [1] TRUE

# It is the derivative of the log-density, so a central difference of the
# log-density in y reproduces it.
eps <- 1e-6
all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
           distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
          distrib_grad_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE

# It vanishes at the mode, which exists because the shape exceeds one.
mode <- 2 * ((1.5 - 1) / 1.5)^(1 / 1.5)
c(mode = mode, deriv = distrib_grad_y(d, mode, th))
#>      mode     deriv 
#> 0.9614997 0.0000000 
```
