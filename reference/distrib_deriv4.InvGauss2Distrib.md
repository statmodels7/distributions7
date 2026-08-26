# Inverse Gaussian Fourth-Order Derivatives in Mean and Shape

Computes the five distinct fourth derivatives of the inverse Gaussian
log-density with respect to \\\mu\\ and \\\lambda\\, in closed form. As
at third order the log-density is linear in \\\lambda\\ apart from
\\\tfrac12\log\lambda\\, so every component naming the shape twice or
more alongside the mean is exactly zero, and the pure shape component is
\\-3/\lambda^4\\.

`expected` is passed straight to the kernel, which substitutes
\\\mathbb{E}\[Y\] = \mu\\; both routes are closed form, so `approx` and
`nsim` are ignored.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_lambda`,
`mu_mu_lambda_lambda`, `mu_lambda_lambda_lambda` and
`lambda_lambda_lambda_lambda`, each of length
`max(length(y), length(mu), length(lambda))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. \\\lambda\\ names this
family's shape parameter.

## See also

[`distrib_deriv3.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGauss2Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss2Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- invgauss2_distrib()
y <- c(1, 2, 3)
th <- list(mu = 2, lambda = 3)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"                 "mu_mu_mu_lambda"            
#> [3] "mu_mu_lambda_lambda"         "mu_lambda_lambda_lambda"    
#> [5] "lambda_lambda_lambda_lambda"

# Two components are exactly zero and the pure shape one is -3/lambda^4.
c(unique(d4$mu_mu_lambda_lambda),
  unique(d4$mu_lambda_lambda_lambda),
  unique(d4$lambda_lambda_lambda_lambda))
#> [1]  0.00000000  0.00000000 -0.03703704
-3 / 3^4
#> [1] -0.03703704

# Expected values: y is replaced by mu.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] -3.375
#> 
#> $mu_mu_mu_lambda
#> [1] 0.375
#> 
#> $mu_mu_lambda_lambda
#> [1] 0
#> 
#> $mu_lambda_lambda_lambda
#> [1] 0
#> 
#> $lambda_lambda_lambda_lambda
#> [1] -0.03703704
#> 

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 2 + eps, lambda = 3))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 2 - eps, lambda = 3))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
#> [1] TRUE
```
