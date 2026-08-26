# Laplace Fourth-Order Derivatives, Rate Parametrization

Computes the five distinct fourth derivatives of the Laplace log-density
with respect to \\\mu\\ and \\\lambda\\. As at third order, every
component involving \\\mu\\ is zero and the only survivor comes from
\\\log\lambda\\: \$\$\dfrac{\partial^4 \ell}{\partial \lambda^4} =
-\dfrac{6}{\lambda^4}.\$\$

None of these depends on the data, so the observed and expected values
coincide: `expected`, `approx` and `nsim` are accepted and have no
effect.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `lambda`. `mu` is not read.
  `lambda` must be strictly positive.

- expected:

  Ignored: the observed values do not depend on the data, so they are
  already their own expectations. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, for the same reason as `expected`.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_lambda`,
`mu_mu_lambda_lambda`, `mu_lambda_lambda_lambda` and
`lambda_lambda_lambda_lambda`, each of length `length(y)`. The first
four are zero and the last is constant at \\-6/\lambda^4\\.

## Notation

\\\ell^{(i j k l)}\\ is the fourth derivative of the log-density with
respect to parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized
superscripts name derivatives.

## See also

[`distrib_deriv3.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Laplace2Distrib.md)
for the order below;
[`distrib_deriv4.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LaplaceDistrib.md)
for the scale parametrization;
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, lambda = 2)

# Only the pure-rate component survives.
vapply(distrib_deriv4(d, y, th), function(v) v[1], numeric(1))
#>                 mu_mu_mu_mu             mu_mu_mu_lambda 
#>                       0.000                       0.000 
#>         mu_mu_lambda_lambda     mu_lambda_lambda_lambda 
#>                       0.000                       0.000 
#> lambda_lambda_lambda_lambda 
#>                      -0.375 
-6 / 2^4
#> [1] -0.375

# A central difference of the third order reproduces it.
eps <- 1e-5
up <- distrib_deriv3(d, y, list(mu = 0.4, lambda = 2 + eps))
dn <- distrib_deriv3(d, y, list(mu = 0.4, lambda = 2 - eps))
all.equal((up$lambda_lambda_lambda - dn$lambda_lambda_lambda) / (2 * eps),
          distrib_deriv4(d, y, th)$lambda_lambda_lambda_lambda,
          tolerance = 1e-6)
#> [1] TRUE
```
