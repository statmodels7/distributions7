# Beta First Derivative in the Response, Mean and Precision

Computes the first derivative of the beta log-density with respect to
the response, in closed form at the implied shapes \\\alpha = \mu\phi\\
and \\\beta = (1-\mu)\phi\\: \$\$\dfrac{\partial \ell}{\partial y} =
\dfrac{\alpha - 1}{y} - \dfrac{\beta - 1}{1 - y}.\$\$ Where both shapes
exceed one it changes sign at the mode \\(\alpha-1)/(\alpha+\beta-2)\\,
so it is positive below the mode and negative above it. Where a shape
falls below one the density is unbounded at the corresponding endpoint
and the derivative does not change sign there.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. An endpoint makes the
  value infinite unless the corresponding shape is exactly 1.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(phi))`,
one value per observation.

## See also

[`distrib_hess_y.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Beta1Distrib.md)
for the second derivative in the response,
[`distrib_gradient.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta1Distrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)

# Written out at the implied shapes.
a <- 0.4 * 5
b <- 0.6 * 5
all.equal(distrib_grad_y(d, y, th), (a - 1) / y - (b - 1) / (1 - y))
#> [1] TRUE

# Zero at the mode (a - 1)/(a + b - 2), positive below and negative above.
mode <- (a - 1) / (a + b - 2)
c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#>         mode      at_mode 
#> 3.333333e-01 4.440892e-16 

# At mu = 1/2 and phi = 2 both shapes are 1 and the density is flat.
distrib_grad_y(d, y, list(mu = 0.5, phi = 2))
#> [1] 0 0 0

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]  2.50 -2.00 -8.75
```
