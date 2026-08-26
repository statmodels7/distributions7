# Lognormal Third-Order Derivatives

Computes the four distinct third derivatives of the lognormal
log-density with respect to \\\mu\\ and \\\sigma^2\\, in closed form.
With \\r = \log y - \mu\\, \$\$\ell^{(\mu\mu\mu)} = 0, \qquad
\ell^{(\mu\mu\sigma^2)} = \dfrac{1}{\sigma^4}, \qquad
\ell^{(\mu\sigma^2\sigma^2)} = \dfrac{2r}{\sigma^6}, \qquad
\ell^{(\sigma^2\sigma^2\sigma^2)} = \dfrac{3r^2}{\sigma^8} -
\dfrac{1}{\sigma^6}.\$\$ The log-density is quadratic in \\\mu\\, so the
third derivative there is zero. As at the lower orders these are the
Gaussian's read at \\\log y\\.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\r\\ with 0 and \\r^2\\ with \\\sigma^2\\: the component odd
in \\r\\ vanishes and the pure variance one becomes \\2/\sigma^6\\. Both
routes are closed form, so no quadrature is run and `approx` and `nsim`
are ignored.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- y:

  A numeric vector of strictly positive observations. With
  `expected = TRUE` only its length is used.

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

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts name
derivatives.

## See also

[`distrib_hessian.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Lognormal1Distrib.md)
for the order below and
[`distrib_deriv4.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Lognormal1Distrib.md)
for the order above;
[`distrib_deriv3.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian2Distrib.md),
which returns the same numbers at \\\log y\\;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)
d3 <- distrib_deriv3(d, y, th)

# Quadratic in mu, so the third derivative there is 0.
d3$mu_mu_mu
#> [1] 0 0 0

# The other three, written out on the log scale.
r <- log(y) - 0.5
all.equal(d3$mu_mu_sigma2, rep(1 / 0.36^2, 3))
#> [1] TRUE
all.equal(d3$mu_sigma2_sigma2, 2 * r / 0.36^3)
#> [1] TRUE

# Expected values: the component odd in r vanishes.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_sigma2
#> [1] 7.716049
#> 
#> $mu_sigma2_sigma2
#> [1] 0
#> 
#> $sigma2_sigma2_sigma2
#> [1] 42.86694
#> 
2 / 0.36^3
#> [1] 42.86694

# Identical to the Gaussian's at log y.
all.equal(d3, distrib_deriv3(gaussian2_distrib(), log(y), th))
#> [1] TRUE
```
