# Gaussian Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the Gaussian
log-density with respect to \\\mu\\ and \\\sigma\\, in closed form.
Writing \\r = y - \mu\\, \$\$\ell^{(\mu\mu\mu\mu)} =
\ell^{(\mu\mu\mu\sigma)} = 0, \qquad \ell^{(\mu\mu\sigma\sigma)} =
-\dfrac{6}{\sigma^4}, \qquad \ell^{(\mu\sigma\sigma\sigma)} = -\dfrac{24
r}{\sigma^5}, \qquad \ell^{(\sigma\sigma\sigma\sigma)} =
\dfrac{6}{\sigma^4} - \dfrac{60 r^2}{\sigma^6}.\$\$ With
`expected = TRUE` the expectations are returned, obtained by replacing
\\r\\ with 0 and \\r^2\\ with \\\sigma^2\\, which leaves
\\-54/\sigma^4\\ in the last component and zero in the two odd ones.
Both routes are closed form, so `approx` and `nsim` are ignored.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma`,
`mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
`sigma_sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell^{(i j k l)}\\ is the fourth derivative of the log-density with
respect to parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized
superscripts name derivatives.

## See also

[`distrib_deriv3.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian1Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian1Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# Quadratic in mu, so the two components with three or more mu are zero.
c(d4$mu_mu_mu_mu[1], d4$mu_mu_mu_sigma[1])
#> [1] 0 0

# The mixed second-second component is constant at -6/sigma^4.
all.equal(d4$mu_mu_sigma_sigma, rep(-6 / 1.5^4, 3))
#> [1] TRUE

# Expected values: -54/sigma^4 in the pure sigma component.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_mu_sigma
#> [1] 0
#> 
#> $mu_mu_sigma_sigma
#> [1] -1.185185
#> 
#> $mu_sigma_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma_sigma
#> [1] -10.66667
#> 
-54 / 1.5^4
#> [1] -10.66667
```
