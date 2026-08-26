# Gaussian Fourth-Order Derivatives in Mean and Variance

Computes the five distinct fourth derivatives of the Gaussian
log-density with respect to \\\mu\\ and \\v = \sigma^2\\, in closed
form. Writing \\r = y - \mu\\, \$\$\ell^{(\mu\mu\mu\mu)} =
\ell^{(\mu\mu\mu v)} = 0, \qquad \ell^{(\mu\mu vv)} = -\dfrac{2}{v^3},
\qquad \ell^{(\mu vvv)} = -\dfrac{6r}{v^4}, \qquad \ell^{(vvvv)} =
\dfrac{3}{v^4} - \dfrac{12 r^2}{v^5}.\$\$ With `expected = TRUE` the
expectations are returned, obtained by replacing \\r\\ with 0 and
\\r^2\\ with \\v\\, which leaves \\-9/v^4\\ in the last component and
zero in the two odd ones. Both routes are closed form, so `approx` and
`nsim` are ignored.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma2`,
`mu_mu_sigma2_sigma2`, `mu_sigma2_sigma2_sigma2` and
`sigma2_sigma2_sigma2_sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density with
respect to parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized
superscripts name derivatives.

## See also

[`distrib_deriv3.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian2Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian2Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- gaussian2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, sigma2 = 4)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"                 "mu_mu_mu_sigma2"            
#> [3] "mu_mu_sigma2_sigma2"         "mu_sigma2_sigma2_sigma2"    
#> [5] "sigma2_sigma2_sigma2_sigma2"

# Quadratic in mu, so the two components with three or more mu are zero.
c(d4$mu_mu_mu_mu[1], d4$mu_mu_mu_sigma2[1])
#> [1] 0 0

# The mixed second-second component is constant at -2/v^3.
all.equal(d4$mu_mu_sigma2_sigma2, rep(-2 / 4^3, 3))
#> [1] TRUE

# Expected values: -9/v^4 in the pure variance component.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_mu_sigma2
#> [1] 0
#> 
#> $mu_mu_sigma2_sigma2
#> [1] -0.03125
#> 
#> $mu_sigma2_sigma2_sigma2
#> [1] 0
#> 
#> $sigma2_sigma2_sigma2_sigma2
#> [1] -0.03515625
#> 
-9 / 4^4
#> [1] -0.03515625
```
