# Poisson Observed Hessian

Computes the second derivative of the Poisson log-mass with respect to
the mean, one value per observation, in closed form:
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}.\$\$ It
depends on the data only through \\y\\, and is exactly zero at an
observed count of zero, because at \\y = 0\\ the log-mass \\-\mu\\ is
linear in \\\mu\\. The term is formed by dividing \\y\\ by \\\mu\\ once
per power, and a zero count contributes zero, so the value at \\y = 0\\
stays zero however small \\\mu\\ is.

On the **link** scale with the default logarithm the chain rule gives
\\\partial^2\ell/\partial\eta^2 = -\mu\\, which carries no data; this is
the defining property of a canonical link, and in exact arithmetic it
makes the observed and the expected information the same matrix there
(see
[`distrib_expected_hessian.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md)).
The generic forms the link-scale value as \\-y + (y - \mu)\\ from the
derivatives on the parameter scale, and for a count much larger than the
mean the two terms cancel to a value much smaller than either. The
relative error of the observed value is then of the order of machine
epsilon times \\y/\mu\\ (about \\4 \times 10^{-3}\\ at \\y = 1000\\,
\\\mu = 10^{-10}\\), while the expected value is computed directly and
stays exact.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- y:

  A numeric vector of counts.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu`, of length
`max(length(y), length(mu))`.

## See also

[`distrib_gradient.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PoissonDistrib.md)
for the score,
[`distrib_expected_hessian.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md)
for the expectation of this quantity, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()
y <- c(0, 2, 7)
th <- list(mu = 3)

# The closed form, written out; zero at an observed zero.
all.equal(distrib_hessian(d, y, th)$mu_mu, -y / 3^2)
#> [1] TRUE

# On the canonical log link the curvature carries no data: one value,
# repeated, equal to -mu.
distrib_hessian(d, y, th, scale = "link")$mu_mu
#> [1] -3 -3 -3

# A central difference of the score reproduces the parameter-scale value.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 3 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 3 - eps))$mu
all.equal((up - dn) / (2 * eps), distrib_hessian(d, y, th)$mu_mu,
          tolerance = 1e-6)
#> [1] TRUE
```
