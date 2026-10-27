# Gaussian Third-Order Derivatives

Computes the four distinct third derivatives of the Gaussian log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form. Writing \\r =
y - \mu\\, \$\$\ell^{(\mu\mu\mu)} = 0, \qquad \ell^{(\mu\mu\sigma)} =
\dfrac{2}{\sigma^3}, \qquad \ell^{(\mu\sigma\sigma)} =
\dfrac{6r}{\sigma^4}, \qquad \ell^{(\sigma\sigma\sigma)} =
-\dfrac{2}{\sigma^3} + \dfrac{12 r^2}{\sigma^5}.\$\$ With
`expected = TRUE` the expectations are returned, obtained by replacing
\\r\\ with 0 and \\r^2\\ with \\\sigma^2\\:
\\\ell^{(\mu\sigma\sigma)}\\, the one component odd in \\r\\, has
expectation 0, the two constant components keep their values and
\\\ell^{(\sigma\sigma\sigma)}\\ becomes \\10/\sigma^3\\. Both routes are
closed form and the method runs no quadrature.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length is used.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned in place of the observed values. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  a fixed number of observations, set per kernel, the kernel runs
  sequentially. Defaults to `1L`.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(i j k)}\\ is the third derivative of the log-density with
respect to parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts
name derivatives; a subscript on \\\ell\\ never does.

## See also

[`distrib_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian1Distrib.md)
for the order below and
[`distrib_deriv4.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian1Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic and for the numerical route that a family without a
closed form uses.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
d3 <- distrib_deriv3(d, y, th)

# The log-density is quadratic in mu, so the third derivative there is 0.
d3$mu_mu_mu
#> [1] 0 0 0

# Two of the other three components are checked against their closed forms.
all.equal(d3$mu_mu_sigma, rep(2 / 1.5^3, 3))
#> [1] TRUE
all.equal(d3$mu_sigma_sigma, 6 * (y - 0.4) / 1.5^4)
#> [1] TRUE

# Expected values: the components odd in the residual vanish.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_sigma
#> [1] 0.5925926
#> 
#> $mu_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma
#> [1] 2.962963
#> 

# A central difference of the Hessian reproduces the same component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_sigma, tolerance = 1e-6)
#> [1] TRUE
```
