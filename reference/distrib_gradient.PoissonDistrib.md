# Poisson Score

Computes the first derivative of the Poisson log-mass with respect to
the mean, one value per observation, in closed form: \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y}{\mu} - 1 = \dfrac{y - \mu}{\mu}.\$\$ The
log-mass is \\y\log\mu - \mu - \log y!\\, so the score is the residual
divided by the mean and its sum vanishes exactly at \\\hat\mu = \bar
y\\.

On the **link** scale with the default logarithm the generic's chain
rule gives \\\partial\ell/\partial\eta = y - \mu\\, the raw residual:
the log is the canonical link of this family and the score there is the
sufficient statistic minus its mean.

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

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of one numeric vector, `mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell\\ is the log-mass of one observation and \\\mu \> 0\\ the mean,
which is also the variance. \\\eta = \log\mu\\ is the parameter on the
link scale. A **canonical** link is the one for which the log-mass is
linear in the sufficient statistic times \\\eta\\.

## See also

[`distrib_hessian.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PoissonDistrib.md)
for the second derivative,
[`distrib_expected_hessian.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md)
for the information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()
y <- c(0, 2, 7)
th <- list(mu = 3)

# The closed form, written out.
all.equal(distrib_gradient(d, y, th)$mu, (y - 3) / 3)
#> [1] TRUE

# On the canonical log link the score is the raw residual.
distrib_gradient(d, y, th, scale = "link")$mu
#> [1] -3 -1  4
y - 3
#> [1] -3 -1  4

# The summed score vanishes at the sample mean, which is the estimate.
set.seed(4)
z <- distrib_rng(d, 1000, list(mu = 4.2))
sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
#> [1] 7.651518e-14
```
