# Laplace Score, Rate Parametrization

Computes the first derivatives of the Laplace log-density with respect
to \\\mu\\ and \\\lambda\\, one value per observation, in closed form.
Writing \\r = y - \mu\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\lambda\\\mathrm{sign}(r), \qquad \dfrac{\partial \ell}{\partial
\lambda} = \dfrac{1}{\lambda} - \|r\|.\$\$

Both are simpler than their scale-parametrization counterparts: the rate
is the natural parameter of the exponential family in \\\|r\|\\, so its
score is the difference between \\1/\lambda\\ and the sufficient
statistic. The score in \\\mu\\ again carries only the sign of the
residual, so the estimate of \\\mu\\ is the sample median.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `lambda` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `lambda`, each of length
`max(length(y), length(mu), length(lambda))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\lambda \> 0\\ the rate, with variance \\2/\lambda^2\\. \\r = y-\mu\\
is the residual. Here \\\lambda\\ is a rate; the same letter names a
penalty parameter above, and the two meet in the lasso, which is this
family with \\\mu\\ held at zero.

## See also

[`distrib_gradient.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LaplaceDistrib.md)
for the scale parametrization,
[`distrib_expected_hessian.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Laplace2Distrib.md)
for the information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, lambda = 2)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
r <- y - 0.4
all.equal(g$mu, 2 * sign(r))
#> [1] TRUE
all.equal(g$lambda, 1 / 2 - abs(r))
#> [1] TRUE

# The rate's score vanishes summed at 1/mean|r|, its estimate.
set.seed(12)
z <- distrib_rng(d, 2000, list(mu = 3, lambda = 0.5))
lam_hat <- 1 / mean(abs(z - median(z)))
sum(distrib_gradient(d, z, list(mu = median(z), lambda = lam_hat))$lambda)
#> [1] -1.709743e-14
```
