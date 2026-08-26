# Generalized Pareto Expected Information

Returns Smith's (1985) closed form, valid for \\\xi \> -1/2\\:
\$\$E\\\left\[\dfrac{\partial^2\ell}{\partial\sigma^2}\right\] =
\dfrac{-1}{(1+2\xi)\sigma^2}, \qquad
E\\\left\[\dfrac{\partial^2\ell}{\partial\sigma\partial\xi}\right\] =
\dfrac{-1}{(1+2\xi)\sigma(1+\xi)}, \qquad
E\\\left\[\dfrac{\partial^2\ell}{\partial\xi^2}\right\] =
\dfrac{-2}{(1+2\xi)(1+\xi)}.\$\$ Every entry is free of the data, so
`approx` and `nsim` are ignored and
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
answers `TRUE`.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector. Its values do not enter the result, which is an
  expectation; only its length does, through recycling.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  Ignored: the expectation is closed form. Accepted so that the
  signature matches the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`.

## Value

A named list of three numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `sigma_sigma`, `xi_xi`, `sigma_xi`. Every entry is `NA` where
`xi <= -0.5`.

## Where it stops existing

At \\\xi \le -1/2\\ the information does **not exist** and every
component is `NA`. The condition is exactly that the integrand be
integrable: written on the probability scale, the second derivative
grows like \\(1-u)^{-2\|\xi\|}\\ near the upper endpoint, which is
integrable if and only if \\\|\xi\| \< 1/2\\. Below that point the
classical asymptotics of the maximum likelihood estimator do not hold
either, so the `NA` reports a quantity that does not exist.

Approaching the boundary the information diverges: measured at \\\sigma
= 1.5\\, the \\\xi\\ component is \\-2.00\\ at \\\xi = 0\\, \\-7.14\\ at
\\-0.3\\ and \\-196.1\\ at \\-0.49\\.

## Checking it

A Monte Carlo average of the observed Hessian is the weaker reference
here, and disagreed with this formula by 9 per cent at \\\xi = -0.3\\
while the formula was right: the second derivative blows up at the upper
endpoint, so the sample mean converges slowly. Integrating on the
**probability** scale instead, \\E\[h\] = \int_0^1 h(Q(u))\\du\\, turns
the endpoint into an ordinary point and agrees to \\10^{-11}\\.

## References

Smith, R. L. (1985). Maximum likelihood estimation in a class of
nonregular cases. *Biometrika* 72, 67-90.

## See also

[`distrib_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GPDDistrib.md)
for the observed curvature,
[`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md)
for the endpoint that causes the condition, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
th <- list(sigma = 1.5, xi = 0.3)
e <- distrib_expected_hessian(d, 0, th)

# Smith's formula written out.
all.equal(unlist(e[c("sigma_sigma", "sigma_xi", "xi_xi")], use.names = FALSE),
          c(-1 / (1.6 * 1.5^2), -1 / (1.6 * 1.5 * 1.3), -2 / (1.6 * 1.3)))
#> [1] TRUE

# It diverges as the shape approaches -1/2, and does not exist at or
# below it.
t(vapply(c(-0.7, -0.5, -0.49, -0.3, 0),
         function(x) {
           u <- distrib_expected_hessian(d, 0, list(sigma = 1.5, xi = x))
           c(xi = x, sigma_sigma = u$sigma_sigma, xi_xi = u$xi_xi)
         }, numeric(3)))
#>         xi sigma_sigma       xi_xi
#> [1,] -0.70          NA          NA
#> [2,] -0.50          NA          NA
#> [3,] -0.49 -22.2222222 -196.078431
#> [4,] -0.30  -1.1111111   -7.142857
#> [5,]  0.00  -0.4444444   -2.000000

# The strategy argument is ignored, the expectation being closed form.
identical(distrib_expected_hessian(d, 0, th),
          distrib_expected_hessian(d, 0, th, approx = "mc", nsim = 10))
#> [1] TRUE
```
