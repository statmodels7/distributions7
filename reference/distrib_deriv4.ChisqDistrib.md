# Chi-Squared Fourth-Order Derivative

Computes the single fourth derivative of the chi-squared log-density
with respect to \\\mu\\, in closed form: \$\$\ell^{(\mu\mu\mu\mu)} =
-\dfrac{\psi'''(\mu/2)}{16},\$\$ with \\\psi'''\\ the third derivative
of the digamma function, one case of \\\ell^{(k)} =
-\psi^{(k-2)}(\mu/2)/2^{k}\\.

The value is free of the response, so `expected`, `approx` and `nsim`
are all without effect and the result is identical to the bit whichever
is asked for.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`, recycled if of length 1. It must be strictly
  positive.

- expected:

  Logical of length 1, and without effect here, the observed and
  expected fourth derivatives being the same number. Defaults to
  `FALSE`.

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

A named list with one numeric vector, `mu_mu_mu_mu`, of length
`max(length(y), length(mu))` and constant within itself when the
parameter is.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. \\\psi\\ is the digamma
function and \\\psi^{(m)}\\ its \\m\\th derivative.

## See also

[`distrib_deriv3.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ChisqDistrib.md)
for the order below and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
y <- c(1, 4, 9)
th <- list(mu = 4)

# A constant, written out with the third derivative of the digamma.
unique(distrib_deriv4(d, y, th)$mu_mu_mu_mu)
#> [1] -0.03087121
-psigamma(2, deriv = 3) / 16
#> [1] -0.03087121

# Free of the response, so asking for the expectation changes nothing.
identical(distrib_deriv4(d, y, th),
          distrib_deriv4(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the third order reproduces it.
eps <- 1e-5
up <- distrib_deriv3(d, y, list(mu = 4 + eps))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 4 - eps))$mu_mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE
```
