# Beta First Derivative in the Response, the Shapes

Computes the first derivative of the beta log-density with respect to
the response, in closed form: \$\$\dfrac{\partial \ell}{\partial y} =
\dfrac{\alpha - 1}{y} - \dfrac{\beta - 1}{1 - y}.\$\$ Where both shapes
exceed one it changes sign at the mode \\(\alpha-1)/(\alpha+\beta-2)\\,
so it is positive below the mode and negative above it. Where a shape
falls below one the density is unbounded at the corresponding endpoint
and the derivative does not change sign there. At \\\alpha = \beta = 1\\
it is identically zero, the density being the uniform.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. An endpoint makes the
  value infinite unless the corresponding shape is exactly 1.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length
`max(length(y), length(alpha), length(beta))`, one value per
observation.

## See also

[`distrib_hess_y.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Beta2Distrib.md)
for the second derivative in the response,
[`distrib_gradient.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta2Distrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)

all.equal(distrib_grad_y(d, y, th), (2 - 1) / y - (5 - 1) / (1 - y))
#> [1] TRUE

# Zero at the mode (alpha - 1)/(alpha + beta - 2).
mode <- (2 - 1) / (2 + 5 - 2)
c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#>    mode at_mode 
#>     0.2     0.0 

# Identically zero at alpha = beta = 1, where the density is the uniform.
distrib_grad_y(d, y, list(alpha = 1, beta = 1))
#> [1] 0 0 0

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]   5.555556  -2.380952 -11.904762
```
