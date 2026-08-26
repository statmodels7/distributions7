# Gaussian Score

Computes the first derivatives of the Gaussian log-density with respect
to \\\mu\\ and \\\sigma\\, one value per observation, in closed form:
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2},
\qquad \dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 -
\sigma^2}{\sigma^3}.\$\$ The arithmetic runs in a compiled kernel
decomposed over the elements of the output, so the result does not
depend on the thread count.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning, giving \\\partial \ell / \partial
\eta_j = h_j'(\eta_j)\\ \partial \ell / \partial \theta_j\\. This method
always returns the parameter scale; the transformation happens in the
generic.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

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

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\sigma \> 0\\ the standard deviation. \\\eta_j\\ is the coordinate of
parameter \\j\\ on the unconstrained scale of its link, and \\h_j' =
\partial \theta_j / \partial \eta_j\\ the chain-rule factor onto it.

## See also

[`distrib_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian1Distrib.md)
for their expectation,
[`distrib_grad_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gaussian1Distrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
all.equal(g$mu, (y - 0.4) / 1.5^2)
#> [1] TRUE
all.equal(g$sigma, ((y - 0.4)^2 - 1.5^2) / 1.5^3)
#> [1] TRUE

# The summed score vanishes at the maximum likelihood estimate.
set.seed(1)
z <- rnorm(200, mean = 3, sd = 2)
mle <- list(mu = mean(z), sigma = sqrt(mean((z - mean(z))^2)))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>            mu         sigma 
#> -5.695531e-15  7.188694e-15 

# On the link scale the sigma component is multiplied by h' = sigma,
# the derivative of the inverse log link; mu rides the identity and is
# unchanged.
distrib_gradient(d, y, th, scale = "link")$sigma / g$sigma
#> [1] 1.5 1.5 1.5
```
