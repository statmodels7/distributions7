# Gaussian Score in Mean and Precision

Computes the first derivatives of the Gaussian log-density with respect
to \\\mu\\ and \\\tau\\, one value per observation, in closed form. With
\\r = y - \mu\\, \$\$\dfrac{\partial \ell}{\partial \mu} = \tau r,
\qquad \dfrac{\partial \ell}{\partial \tau} = \dfrac{1}{2\tau} -
\dfrac{r^2}{2}.\$\$ The precision component is linear in \\r^2\\. That
linearity is why every third and fourth derivative of this
parametrization is free of the response. The arithmetic runs in a
compiled kernel decomposed over the elements of the output, so the
result does not depend on the thread count.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning, giving \\\partial \ell / \partial
\eta_j = h_j'(\eta_j)\\ \partial \ell / \partial \theta_j\\. This method
always returns the parameter scale; the transformation happens in the
generic.

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

A named list of two numeric vectors, `mu` and `tau`, each of length
`max(length(y), length(mu), length(tau))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\tau = 1/\sigma^2 \> 0\\ the precision. \\\eta_j\\ is the coordinate
of parameter \\j\\ on the unconstrained scale of its link, and \\h_j' =
\partial \theta_j / \partial \eta_j\\ the chain-rule factor onto it.

## See also

[`distrib_hessian.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian3Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian3Distrib.md)
for their expectation,
[`distrib_gradient.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian1Distrib.md)
for the same score in the standard deviation, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, tau = 0.25)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
r <- y - 1
all.equal(g$mu, 0.25 * r)
#> [1] TRUE
all.equal(g$tau, 1 / (2 * 0.25) - r^2 / 2)
#> [1] TRUE

# The summed score vanishes at the maximum likelihood estimate, where the
# precision is the reciprocal of the variance with divisor n.
set.seed(1)
z <- rnorm(200, mean = 3, sd = 2)
mle <- list(mu = mean(z), tau = 1 / mean((z - mean(z))^2))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>            mu           tau 
#> -5.859029e-15  5.218048e-14 

# log(tau) is minus twice log(sigma), so the link-scale score here is minus
# half gaussian1's.
distrib_gradient(d, y, th, scale = "link")$tau /
  distrib_gradient(gaussian1_distrib(), y, list(mu = 1, sigma = 2),
                   scale = "link")$sigma
#> [1] -0.5 -0.5 -0.5
```
