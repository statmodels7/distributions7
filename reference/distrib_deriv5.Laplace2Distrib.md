# Laplace, Rate Parametrization Fifth-Order Derivatives

Computes the six distinct fifth derivatives of the Laplace log-density
in the rate parametrization, \\\ell = \log(\lambda/2) - \lambda\|y -
\mu\|\\, in closed form from a compiled kernel. Away from the kink at
\\y = \mu\\ only the pure \\\lambda\\ component is nonzero,
\$\$\ell^{(\lambda\lambda\lambda\lambda\lambda)} =
\dfrac{24}{\lambda^5}.\$\$ With `scale = "link"` the components are
composed with the link derivatives by the routine that every family with
a fifth-order kernel uses.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`.

- scale:

  One of `"parameter"` (the default) or `"link"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of six numeric vectors, `mu_mu_mu_mu_mu`,
`mu_mu_mu_mu_lambda`, `mu_mu_mu_lambda_lambda`,
`mu_mu_lambda_lambda_lambda`, `mu_lambda_lambda_lambda_lambda` and
`lambda_lambda_lambda_lambda_lambda`, each of length
`max(length(y), length(mu), length(lambda))`.

## See also

[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
d5 <- distrib_deriv5(d, c(-1.2, 0.3, 2.5), list(mu = 0.4, lambda = 2))
all.equal(d5$lambda_lambda_lambda_lambda_lambda, rep(24 / 2^5, 3))
#> [1] TRUE
```
