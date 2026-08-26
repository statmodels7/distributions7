# Student t Second Derivative in the Response

Computes \\\partial^2 \ell / \partial y^2\\, the second derivative of
the Student t log-density with respect to the response, in closed form.
With \\r = y - \mu\\ and \\D = \nu\sigma^2 + r^2\\,
\$\$\dfrac{\partial^2 \ell}{\partial y^2} = \dfrac{(\nu+1)\left(r^2 -
\nu\sigma^2\right)}{D^2}.\$\$ The family is a location family, so this
equals the pure-location entry of
[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md),
with no sign change; two derivatives in \\\mu\\ carry two factors of
\\-1\\. It is negative near the mode and **positive** beyond \\\|r\| =
\sigma\sqrt{\nu}\\, so the log-density is not concave in the response,
unlike a Gaussian's.

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

[`distrib_grad_y.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.StudentT1Distrib.md)
for the first derivative in the response,
[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md)
for the second derivatives in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)

# The closed form, written out.
r <- y - 0.4; D <- 5 * 1.2^2 + r^2
all.equal(distrib_hess_y(d, y, th), 6 * (r^2 - 5 * 1.2^2) / D^2)
#> [1] TRUE

# A location family, so this is the pure-location entry of the Hessian.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# A central difference of the first derivative reproduces it.
eps <- 1e-5
all.equal((distrib_grad_y(d, y + eps, th) -
           distrib_grad_y(d, y - eps, th)) / (2 * eps),
          distrib_hess_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE

# Positive beyond |r| = sigma * sqrt(nu) = 2.68, so the log-density is
# convex in the response out there.
distrib_hess_y(d, 0.4 + c(1, 2, 4, 8), th)
#> [1] -0.55324212 -0.15306122  0.09809750  0.06722636
```
