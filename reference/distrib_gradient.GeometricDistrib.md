# Geometric Score

Computes the first derivative of the geometric log-mass with respect to
the mean, one value per observation, in closed form: \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y - \mu}{\mu(1+\mu)}.\$\$ The residual is
divided by the variance \\\mu(1+\mu)\\, and the sum vanishes exactly at
\\\hat\mu = \bar y\\.

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

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

A named list of one numeric vector, `mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell\\ is the log-mass of one observation and \\\mu \> 0\\ the mean,
with variance \\\mu(1+\mu)\\. The success probability is \\p =
1/(1+\mu)\\.

## See also

[`distrib_hessian.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GeometricDistrib.md)
for the second derivative,
[`distrib_expected_hessian.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GeometricDistrib.md)
for the information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- geometric_distrib()
y <- c(0, 2, 7)
th <- list(mu = 3)

# The closed form, written out: the residual over the variance.
all.equal(distrib_gradient(d, y, th)$mu, (y - 3) / (3 * (1 + 3)))
#> [1] TRUE

# The summed score vanishes at the sample mean, which is the estimate.
set.seed(5)
z <- distrib_rng(d, 2000, list(mu = 2.5))
sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
#> [1] 2.519512e-14
```
