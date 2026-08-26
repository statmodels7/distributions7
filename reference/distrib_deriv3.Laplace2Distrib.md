# Laplace Third-Order Derivatives, Rate Parametrization

Computes the four distinct third derivatives of the Laplace log-density
with respect to \\\mu\\ and \\\lambda\\. In this parametrization the
log-density is \\\log(\lambda/2) - \lambda\|y-\mu\|\\, which is **linear
in \\\lambda\\** apart from \\\log\lambda\\, and piecewise linear in
\\\mu\\. Every third derivative involving \\\mu\\ at all is therefore
zero, and the only surviving component comes from \\\log\lambda\\:
\$\$\dfrac{\partial^3 \ell}{\partial \lambda^3} =
\dfrac{2}{\lambda^3},\$\$ with `mu_mu_mu`, `mu_mu_lambda` and
`mu_lambda_lambda` all zero.

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

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. `mu` is not read. `lambda` must
  be strictly positive.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_lambda`,
`mu_lambda_lambda` and `lambda_lambda_lambda`, each of length
`length(y)`. The first three are zero and the last is constant at
\\2/\lambda^3\\.

## Notation

\\\ell^{(i j k)}\\ is the third derivative of the log-density with
respect to parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts
name derivatives.

## See also

[`distrib_deriv3.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LaplaceDistrib.md)
for the scale parametrization, where the picture is less sparse;
[`distrib_deriv4.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Laplace2Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, lambda = 2)
d3 <- distrib_deriv3(d, y, th)

# Only the pure-rate component survives.
vapply(d3, function(v) v[1], numeric(1))
#>             mu_mu_mu         mu_mu_lambda     mu_lambda_lambda 
#>                 0.00                 0.00                 0.00 
#> lambda_lambda_lambda 
#>                 0.25 
2 / 2^3
#> [1] 0.25

# expected = TRUE changes nothing, the values carrying no data.
identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#> [1] TRUE
```
