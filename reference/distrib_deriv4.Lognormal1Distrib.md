# Lognormal Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the lognormal
log-density with respect to \\\mu\\ and \\\sigma^2\\, in closed form.
With \\r = \log y - \mu\\, \$\$\ell^{(\mu\mu\mu\mu)} =
\ell^{(\mu\mu\mu\sigma^2)} = 0, \qquad \ell^{(\mu\mu\sigma^2\sigma^2)} =
-\dfrac{2}{\sigma^6}, \qquad \ell^{(\mu\sigma^2\sigma^2\sigma^2)} =
-\dfrac{6r}{\sigma^8}, \qquad \ell^{(\sigma^{2\cdot 4})} =
\dfrac{3}{\sigma^8} - \dfrac{12 r^2}{\sigma^{10}}.\$\$ As at the lower
orders these are the Gaussian's read at \\\log y\\.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\r\\ with 0 and \\r^2\\ with \\\sigma^2\\, which leaves
\\-9/\sigma^8\\ in the last component and zero in the two odd ones. Both
routes are closed form, so `approx` and `nsim` are ignored.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma2`,
`mu_mu_sigma2_sigma2`, `mu_sigma2_sigma2_sigma2` and
`sigma2_sigma2_sigma2_sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized superscripts
name derivatives.

## See also

[`distrib_deriv3.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Lognormal1Distrib.md)
for the order below,
[`distrib_deriv4.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian2Distrib.md),
which returns the same numbers at \\\log y\\, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"                 "mu_mu_mu_sigma2"            
#> [3] "mu_mu_sigma2_sigma2"         "mu_sigma2_sigma2_sigma2"    
#> [5] "sigma2_sigma2_sigma2_sigma2"

# Quadratic in mu, so the two components with three or more mu are zero.
c(d4$mu_mu_mu_mu[1], d4$mu_mu_mu_sigma2[1])
#> [1] 0 0

# The mixed second-second component is constant at -2/sigma2^3.
all.equal(d4$mu_mu_sigma2_sigma2, rep(-2 / 0.36^3, 3))
#> [1] TRUE

# Expected values: -9/sigma2^4 in the pure variance component.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_mu_sigma2
#> [1] 0
#> 
#> $mu_mu_sigma2_sigma2
#> [1] -42.86694
#> 
#> $mu_sigma2_sigma2_sigma2
#> [1] 0
#> 
#> $sigma2_sigma2_sigma2_sigma2
#> [1] -535.8368
#> 
-9 / 0.36^4
#> [1] -535.8368
```
