# Student t First Derivative in the Response

Computes \\\partial \ell / \partial y\\, the derivative of the Student t
log-density with respect to the response, in closed form. With \\r = y -
\mu\\ and \\D = \nu\sigma^2 + r^2\\, \$\$\dfrac{\partial \ell}{\partial
y} = -\dfrac{(\nu+1)r}{D}.\$\$ The family is a location family in
\\\mu\\, so this is exactly the negative of the location score
[`distrib_gradient.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md)`$mu`,
and it redescends with it: it is largest in size at \\\|r\| =
\sigma\sqrt{\nu}\\ and falls back towards zero beyond that. This
quantity is what a quantile residual's delta-method standard error and a
change of variable in the response both need.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length
`max(length(y), length(mu), length(sigma), length(nu))`, one value per
observation.

## See also

[`distrib_hess_y.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.StudentT1Distrib.md)
for the second derivative in the response,
[`distrib_cross_y.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.StudentT1Distrib.md)
for the mixed derivative in the response and the parameters,
[`distrib_gradient.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md)
for the derivatives in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)

# The closed form, written out.
r <- y - 0.4
all.equal(distrib_grad_y(d, y, th), -6 * r / (5 * 1.2^2 + r^2))
#> [1] TRUE

# A location family, so this is minus the location score.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# It is the derivative of the log-density, so a central difference of the
# log-density in y reproduces it.
eps <- 1e-6
all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
           distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
          distrib_grad_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE
```
