# Chi-Squared First Derivative in the Response

Computes the first derivative of the chi-squared log-density with
respect to the response, in closed form: \$\$\dfrac{\partial
\ell}{\partial y} = \dfrac{\mu/2 - 1}{y} - \dfrac{1}{2}.\$\$ It changes
sign at \\y = \mu - 2\\, which is the mode of the density for \\\mu \>
2\\. At \\\mu = 2\\ the first term drops out and the derivative is the
constant \\-1/2\\ of an exponential with mean 2; below \\\mu = 2\\ the
density has no interior mode and the derivative is negative throughout.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of strictly positive observations. At `y = 0` the
  value is infinite unless \\\mu\\ is exactly 2.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`, recycled if of length 1. It must be strictly
  positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu))`, one value per
observation.

## See also

[`distrib_hess_y.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ChisqDistrib.md)
for the second derivative in the response,
[`distrib_gradient.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ChisqDistrib.md)
for the derivative in the parameter, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
y <- c(1, 4, 9)
th <- list(mu = 4)

# Written out.
all.equal(distrib_grad_y(d, y, th), (4 / 2 - 1) / y - 0.5)
#> [1] TRUE

# Zero at the mode mu - 2, positive below it and negative above.
c(mode = 4 - 2, at_mode = distrib_grad_y(d, 2, th))
#>    mode at_mode 
#>       2       0 
distrib_grad_y(d, c(1, 9), th)
#> [1]  0.5000000 -0.3888889

# At mu = 2 the family is exponential and the derivative is constant.
distrib_grad_y(d, y, list(mu = 2))
#> [1] -0.5 -0.5 -0.5

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]  0.5000000 -0.2500000 -0.3888889
```
