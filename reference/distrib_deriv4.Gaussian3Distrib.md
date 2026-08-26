# Gaussian Fourth-Order Derivatives in Mean and Precision

Computes the five distinct fourth derivatives of the Gaussian
log-density with respect to \\\mu\\ and \\\tau\\, in closed form. Four
of the five are exactly zero and the fifth is a constant:
\$\$\ell^{(\tau\tau\tau\tau)} = -\dfrac{3}{\tau^4}.\$\$ This is the
flattest of the three parametrizations of the Gaussian at fourth order;
the log-density is quadratic in \\\mu\\ and its \\\tau\\ part is
\\\tfrac{1}{2}\log\tau\\ plus a term linear in \\\tau\\, so every mixed
component past \\\ell^{(\mu\mu\tau)}\\ vanishes. Like the third order,
the values are free of the response, so `expected`, `approx` and `nsim`
are all without effect.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- expected:

  Logical of length 1, and without effect here, the observed and
  expected fourth derivatives being the same numbers. Defaults to
  `FALSE`.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of the length of `y`. `mu` is not read. `tau` must be
  strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_tau`,
`mu_mu_tau_tau`, `mu_tau_tau_tau` and `tau_tau_tau_tau`, each of length
`length(y)`. Only the last is non-zero.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density with
respect to parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized
superscripts name derivatives.

## See also

[`distrib_deriv3.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian3Distrib.md)
for the order below,
[`distrib_deriv4.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian1Distrib.md)
for the same order in the standard deviation, where three of the five
are non-zero, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, tau = 0.25)
d4 <- distrib_deriv4(d, y, th)

# One non-zero component out of five, at -3/tau^4.
lapply(d4, unique)
#> $mu_mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_mu_tau
#> [1] 0
#> 
#> $mu_mu_tau_tau
#> [1] 0
#> 
#> $mu_tau_tau_tau
#> [1] 0
#> 
#> $tau_tau_tau_tau
#> [1] -768
#> 
-3 / 0.25^4
#> [1] -768

# Free of the response, so asking for the expectation changes nothing.
identical(d4, distrib_deriv4(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 1, tau = 0.25 + eps))$tau_tau_tau
dn <- distrib_deriv3(d, y, list(mu = 1, tau = 0.25 - eps))$tau_tau_tau
all.equal((up - dn) / (2 * eps), d4$tau_tau_tau_tau, tolerance = 1e-4)
#> [1] TRUE
```
