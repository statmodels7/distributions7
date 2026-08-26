# Logistic Score

Computes the first derivatives of the logistic log-density with respect
to \\\mu\\ and \\\sigma\\, one value per observation, in closed form.
Writing \\z = (y - \mu)/\sigma\\, \$\$\dfrac{\partial \ell}{\partial
\mu} = \dfrac{1}{\sigma} \tanh\left(\dfrac{z}{2}\right), \qquad
\dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1}{\sigma}\left\[1 - z
\tanh\left(\dfrac{z}{2}\right)\right\].\$\$

The score in \\\mu\\ is bounded by \\1/\sigma\\ and saturates rather
than redescending: a distant observation contributes a fixed amount
instead of an unbounded one, as it would under a Gaussian, or a
vanishing one, as it would under a Cauchy.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

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

A named list of two numeric vectors, `mu` and `sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean,
\\\sigma \> 0\\ the scale and \\z = (y-\mu)/\sigma\\ the standardized
residual. \\\sigma\\ is not the standard deviation, which is
\\\pi\sigma/\sqrt{3}\\.

## See also

[`distrib_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LogisticDistrib.md)
for their expectation,
[`distrib_gradient.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md)
for a redescending score and
[`distrib_gradient.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian1Distrib.md)
for an unbounded one, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
z <- (y - 0.4) / 1.5
all.equal(g$mu, tanh(z / 2) / 1.5)
#> [1] TRUE
all.equal(g$sigma, -(1 - z * tanh(z / 2)) / 1.5)
#> [1] TRUE

# The score in mu saturates at 1/sigma instead of growing.
round(distrib_gradient(d, 0.4 + c(0, 1.5, 3, 6, 12, 60), th)$mu, 4)
#> [1] 0.0000 0.3081 0.5077 0.6427 0.6662 0.6667
1 / 1.5
#> [1] 0.6666667

# The summed score vanishes at the maximum likelihood estimate.
set.seed(9)
zz <- distrib_rng(d, 3000, list(mu = 2, sigma = 1))
fit <- fit_distrib(d, zz)
round(vapply(distrib_gradient(d, zz, as.list(coef(fit))), sum, numeric(1)), 8)
#>          mu       sigma 
#>  0.00025606 -0.00002239 
```
