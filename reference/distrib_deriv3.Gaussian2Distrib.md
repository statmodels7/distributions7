# Gaussian Third-Order Derivatives in Mean and Variance

Computes the four distinct third derivatives of the Gaussian log-density
with respect to \\\mu\\ and \\v = \sigma^2\\, in closed form. Writing
\\r = y - \mu\\, \$\$\ell^{(\mu\mu\mu)} = 0, \qquad \ell^{(\mu\mu v)} =
\dfrac{1}{v^2}, \qquad \ell^{(\mu v v)} = \dfrac{2r}{v^3}, \qquad
\ell^{(vvv)} = \dfrac{3r^2}{v^4} - \dfrac{1}{v^3}.\$\$ With
`expected = TRUE` the expectations are returned, obtained by replacing
\\r\\ with 0 and \\r^2\\ with \\v\\: the component odd in \\r\\ vanishes
and \\\ell^{(vvv)}\\ becomes \\2/v^3\\. Both routes are closed form, so
no quadrature is run and `approx` and `nsim` are ignored.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length is used.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma2`,
`mu_sigma2_sigma2` and `sigma2_sigma2_sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density with respect
to parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts name
derivatives; a subscript on \\\ell\\ never does.

## See also

[`distrib_hessian.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian2Distrib.md)
for the order below and
[`distrib_deriv4.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian2Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic and for the numerical route a family without a closed
form takes.

## Examples

``` r
d <- gaussian2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, sigma2 = 4)
d3 <- distrib_deriv3(d, y, th)

# The log-density is quadratic in mu, so the third derivative there is 0.
d3$mu_mu_mu
#> [1] 0 0 0

# The other three, written out.
r <- y - 1
all.equal(d3$mu_mu_sigma2, rep(1 / 4^2, 3))
#> [1] TRUE
all.equal(d3$mu_sigma2_sigma2, 2 * r / 4^3)
#> [1] TRUE
all.equal(d3$sigma2_sigma2_sigma2, 3 * r^2 / 4^4 - 1 / 4^3)
#> [1] TRUE

# Expected values: the component odd in the residual vanishes and the pure
# variance component becomes 2/v^3.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_sigma2
#> [1] 0.0625
#> 
#> $mu_sigma2_sigma2
#> [1] 0
#> 
#> $sigma2_sigma2_sigma2
#> [1] 0.03125
#> 
2 / 4^3
#> [1] 0.03125

# A central difference of the Hessian reproduces the same component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 1, sigma2 = 4 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 1, sigma2 = 4 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_sigma2, tolerance = 1e-6)
#> [1] TRUE
```
