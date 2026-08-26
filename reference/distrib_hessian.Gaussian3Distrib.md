# Gaussian Observed Hessian in Mean and Precision

Computes the three distinct second derivatives of the Gaussian
log-density with respect to \\\mu\\ and \\\tau\\, one value per
observation, in closed form. With \\r = y - \mu\\, \$\$\ell^{(\mu\mu)} =
-\tau, \qquad \ell^{(\mu\tau)} = r, \qquad \ell^{(\tau\tau)} =
-\dfrac{1}{2\tau^2}.\$\$ Two of the three are free of the data, so the
observed and the expected Hessians differ in the mixed entry alone. That
entry sums to \\\sum_i (y_i - \mu)\\, which vanishes at \\\hat\mu = \bar
y\\: the two matrices agree exactly once the mean equation is solved and
differ before that, so Fisher scoring and Newton's method take the same
step at the maximum and different steps on the way there.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `tau` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_tau` and `tau_tau`,
each of length `max(length(y), length(mu), length(tau))`. The three name
the distinct entries of a symmetric \\2 \times 2\\ matrix per
observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density with respect
to parameters \\i\\ and \\j\\. Parenthesized superscripts name
derivatives; a subscript on \\\ell\\ never does.

## See also

[`distrib_gradient.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian3Distrib.md)
for the score,
[`distrib_expected_hessian.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian3Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian3Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, tau = 0.25)
h <- distrib_hessian(d, y, th)

# Two of the three are constants; only the mixed entry carries the data.
h$mu_mu
#> [1] -0.25 -0.25 -0.25
h$tau_tau
#> [1] -8 -8 -8
all.equal(h$mu_tau, y - 1)
#> [1] TRUE

# The mixed entry sums to zero at the estimated mean and not elsewhere.
set.seed(3)
z <- distrib_rng(d, 500, list(mu = 3, tau = 0.25))
c(at_mle = sum(distrib_hessian(d, z, list(mu = mean(z), tau = 0.25))$mu_tau),
  at_2   = sum(distrib_hessian(d, z, list(mu = 2, tau = 0.25))$mu_tau))
#>        at_mle          at_2 
#> -6.705747e-14  5.525412e+02 

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 1, tau = 0.25 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 1, tau = 0.25 - eps))$mu
all.equal((up - dn) / (2 * eps), h$mu_tau, tolerance = 1e-6)
#> [1] TRUE
```
