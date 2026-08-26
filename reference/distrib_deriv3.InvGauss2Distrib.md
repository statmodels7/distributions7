# Inverse Gaussian Third-Order Derivatives in Mean and Shape

Computes the four distinct third derivatives of the inverse Gaussian
log-density with respect to \\\mu\\ and \\\lambda\\, in closed form:
\$\$\ell^{(\mu\mu\mu)} = \dfrac{\lambda(12y - 6\mu)}{\mu^5}, \qquad
\ell^{(\mu\mu\lambda)} = \dfrac{2\mu - 3y}{\mu^4}, \qquad
\ell^{(\mu\lambda\lambda)} = 0, \qquad \ell^{(\lambda\lambda\lambda)} =
\dfrac{1}{\lambda^3}.\$\$ The log-density is linear in \\\lambda\\ apart
from \\\tfrac12\log\lambda\\, so any component naming the shape twice
and the mean once is exactly zero and the pure shape component carries
only the logarithm. `expected` is passed straight to the kernel, which
substitutes \\\mathbb{E}\[Y\] = \mu\\; both routes are closed form, so
`approx` and `nsim` are ignored.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- y:

  A numeric vector of strictly positive observations. With
  `expected = TRUE` only its length is used.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned in place of the observed values. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available for both the observed and the
  expected values.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_lambda`,
`mu_lambda_lambda` and `lambda_lambda_lambda`, each of length
`max(length(y), length(mu), length(lambda))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. \\\lambda\\ names this family's shape
parameter.

## See also

[`distrib_hessian.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss2Distrib.md)
for the order below and
[`distrib_deriv4.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.InvGauss2Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()
y <- c(1, 2, 3)
th <- list(mu = 2, lambda = 3)
d3 <- distrib_deriv3(d, y, th)

# The mixed mu-lambda-lambda component is exactly zero.
d3$mu_lambda_lambda
#> [1] 0 0 0

# The pure shape component is 1/lambda^3, free of the data.
unique(d3$lambda_lambda_lambda)
#> [1] 0.03703704
1 / 3^3
#> [1] 0.03703704

# Expected values: y is replaced by mu.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 1.125
#> 
#> $mu_mu_lambda
#> [1] -0.125
#> 
#> $mu_lambda_lambda
#> [1] 0
#> 
#> $lambda_lambda_lambda
#> [1] 0.03703704
#> 

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 2 + eps, lambda = 3))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 2 - eps, lambda = 3))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
#> [1] TRUE
```
