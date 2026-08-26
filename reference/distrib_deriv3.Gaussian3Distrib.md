# Gaussian Third-Order Derivatives in Mean and Precision

Computes the four distinct third derivatives of the Gaussian log-density
with respect to \\\mu\\ and \\\tau\\, in closed form:
\$\$\ell^{(\mu\mu\mu)} = 0, \qquad \ell^{(\mu\mu\tau)} = -1, \qquad
\ell^{(\mu\tau\tau)} = 0, \qquad \ell^{(\tau\tau\tau)} =
\dfrac{1}{\tau^3}.\$\$ Every one of them is free of the response, so the
observed and the expected values coincide and `expected` selects
nothing: the same kernel runs either way and the two results are
identical to the bit. `approx` and `nsim` are ignored for the same
reason.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- expected:

  Logical of length 1, and without effect here, the observed and
  expected third derivatives being the same numbers. Defaults to
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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_tau`,
`mu_tau_tau` and `tau_tau_tau`, each of length `length(y)`. The names
enumerate the distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density with respect
to parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts name
derivatives; a subscript on \\\ell\\ never does.

## See also

[`distrib_hessian.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian3Distrib.md)
for the order below and
[`distrib_deriv4.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian3Distrib.md)
for the order above;
[`distrib_deriv3.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian1Distrib.md)
for the same order in the standard deviation, where two components do
carry the residual; and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, tau = 0.25)
d3 <- distrib_deriv3(d, y, th)

# Two of the four are non-zero, and both are constants.
lapply(d3, unique)
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_tau
#> [1] -1
#> 
#> $mu_tau_tau
#> [1] 0
#> 
#> $tau_tau_tau
#> [1] 64
#> 
1 / 0.25^3
#> [1] 64

# Free of the response, so asking for the expectation changes nothing.
identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 1, tau = 0.25 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 1, tau = 0.25 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_tau, tolerance = 1e-6)
#> [1] TRUE
```
