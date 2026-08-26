# Inverse Gaussian Score in Mean and Shape

Computes the first derivatives of the inverse Gaussian log-density with
respect to \\\mu\\ and \\\lambda\\, one value per observation, in closed
form: \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\lambda(y-\mu)}{\mu^3}, \qquad \dfrac{\partial \ell}{\partial
\lambda} = \dfrac{1}{2\lambda} - \dfrac{(y-\mu)^2}{2\mu^2 y}.\$\$ The
log-density is **linear in the shape** apart from
\\\tfrac12\log\lambda\\. Every derivative of this family is elementary
for that reason, and every expectation reduces to \\\mathbb{E}\[Y\] =
\mu\\ and \\\mathbb{E}\[(Y-\mu)^2/Y\] = \mu^2/\lambda\\. The arithmetic
runs in a compiled kernel decomposed over the elements of the output, so
the result does not depend on the thread count.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

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

A named list of two numeric vectors, `mu` and `lambda`, each of length
`max(length(y), length(mu), length(lambda))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the mean
and \\\lambda \> 0\\ the shape, with \\\operatorname{Var}(Y) =
\mu^3/\lambda\\. Here \\\lambda\\ names this family's shape parameter
throughout.

## See also

[`distrib_hessian.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss2Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss2Distrib.md)
for their expectation,
[`distrib_gradient.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss1Distrib.md)
for the same score in the dispersion, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()
y <- c(1, 2, 3)
th <- list(mu = 2, lambda = 3)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
all.equal(g$mu, 3 * (y - 2) / 2^3)
#> [1] TRUE
all.equal(g$lambda, 1 / (2 * 3) - (y - 2)^2 / (2 * 2^2 * y))
#> [1] TRUE

# The mean component vanishes at y = mu, whatever the shape.
distrib_gradient(d, 2, list(mu = 2, lambda = c(0.5, 3, 30)))$mu
#> [1] 0 0 0

# Both estimating equations solve in closed form, so the summed score
# vanishes at the sample mean and 1/(mean(1/y) - 1/mean(y)).
set.seed(3)
z <- distrib_rng(d, 2000, th)
mle <- list(mu = mean(z), lambda = 1 / (mean(1 / z) - 1 / mean(z)))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>            mu        lambda 
#> -3.414283e-14  4.635181e-15 
```
