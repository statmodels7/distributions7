# Inverse Gaussian Score in Mean and Dispersion

Computes the first derivatives of the inverse Gaussian log-density with
respect to \\\mu\\ and \\\phi\\, one value per observation, in closed
form: \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y -
\mu}{\phi\mu^3}, \qquad \dfrac{\partial \ell}{\partial \phi} =
\dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2}.\$\$ The mean component
is the score of an inverse Gaussian generalized linear model, the
residual divided by the variance function \\\phi\mu^3\\. The arithmetic
runs in a compiled kernel decomposed over the elements of the output, so
the result does not depend on the thread count.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

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

A named list of two numeric vectors, `mu` and `phi`, each of length
`max(length(y), length(mu), length(phi))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the mean
and \\\phi \> 0\\ the dispersion, with \\\operatorname{Var}(Y) =
\phi\mu^3\\.

## See also

[`distrib_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss1Distrib.md)
for their expectation,
[`distrib_grad_y.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss1Distrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)
g <- distrib_gradient(d, y, th)

# The generalized linear model score: residual over the variance function.
all.equal(g$mu, (y - 1) / (2 * 1^3))
#> [1] TRUE
all.equal(g$phi, ((y - 1)^2 - y * 1^2 * 2) / (2 * y * 2^2 * 1^2))
#> [1] TRUE

# The mean component vanishes at y = mu, whatever the dispersion.
distrib_gradient(d, 1, list(mu = 1, phi = c(0.5, 2, 8)))$mu
#> [1] 0 0 0

# Both estimating equations solve in closed form, so the summed score
# vanishes at the sample mean and mean(1/y) - 1/mean(y).
set.seed(3)
z <- distrib_rng(d, 2000, th)
mle <- list(mu = mean(z), phi = mean(1 / z) - 1 / mean(z))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>            mu           phi 
#> -1.022116e-13  7.635060e-14 
```
