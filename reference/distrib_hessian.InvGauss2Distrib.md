# Inverse Gaussian Observed Hessian in Mean and Shape

Computes the three distinct second derivatives of the inverse Gaussian
log-density with respect to \\\mu\\ and \\\lambda\\, one value per
observation, in closed form: \$\$\ell^{(\mu\mu)} = \dfrac{\lambda(2\mu -
3y)}{\mu^4}, \qquad \ell^{(\mu\lambda)} = \dfrac{y - \mu}{\mu^3}, \qquad
\ell^{(\lambda\lambda)} = -\dfrac{1}{2\lambda^2}.\$\$

The shape parametrization is better behaved than the dispersion one at
second order. The pure shape entry is a negative constant, free of the
data, where the corresponding entry of
[`distrib_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss1Distrib.md)
can be positive; only the curvature in \\\mu\\ can turn positive here,
and it does so wherever \\y \< 2\mu/3\\.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_lambda` and
`lambda_lambda`, in that order, each of length
`max(length(y), length(mu), length(lambda))`. The three name the
distinct entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.
\\\lambda\\ names this family's shape parameter.

## See also

[`distrib_gradient.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss2Distrib.md)
for the score,
[`distrib_expected_hessian.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss2Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGauss2Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()
y <- c(1, 2, 3)
th <- list(mu = 2, lambda = 3)
h <- distrib_hessian(d, y, th)
h
#> $mu_mu
#> [1]  0.1875 -0.3750 -0.9375
#> 
#> $mu_lambda
#> [1] -0.125  0.000  0.125
#> 
#> $lambda_lambda
#> [1] -0.05555556 -0.05555556 -0.05555556
#> 

# The three closed forms, written out.
all.equal(h$mu_mu, 3 * (2 * 2 - 3 * y) / 2^4)
#> [1] TRUE
all.equal(h$mu_lambda, (y - 2) / 2^3)
#> [1] TRUE
all.equal(h$lambda_lambda, rep(-1 / (2 * 3^2), 3))
#> [1] TRUE

# The curvature in mu is positive below 2 mu/3 and negative above.
c(at_1 = distrib_hessian(d, 1, th)$mu_mu,
  at_3 = distrib_hessian(d, 3, th)$mu_mu)
#>    at_1    at_3 
#>  0.1875 -0.9375 

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 2 + eps, lambda = 3))$mu
dn <- distrib_gradient(d, y, list(mu = 2 - eps, lambda = 3))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
