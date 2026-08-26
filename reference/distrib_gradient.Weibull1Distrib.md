# Weibull Score

Computes the first derivatives of the Weibull log-density with respect
to the scale \\\mu\\ and the shape \\\sigma\\, one value per
observation, in closed form. Writing \\z = y/\mu\\ and \\u =
z^{\sigma}\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\sigma}{\mu}(u - 1), \qquad \dfrac{\partial \ell}{\partial
\sigma} = \dfrac{1}{\sigma} + (1 - u)\log z.\$\$ Both components are
polynomials in \\u\\ and \\u \log z\\, so every higher order is
elementary as well.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning, giving \\\partial \ell / \partial
\eta_j = h_j'(\eta_j)\\ \partial \ell / \partial \theta_j\\. This method
always returns the parameter scale; the transformation happens in the
generic.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- y:

  A numeric vector of positive observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the scale
and \\\sigma \> 0\\ the shape. \\\eta_j\\ is the coordinate of parameter
\\j\\ on the unconstrained scale of its link, and \\h_j' = \partial
\theta_j / \partial \eta_j\\ the chain-rule factor onto it.

## See also

[`distrib_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Weibull1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Weibull1Distrib.md)
for their expectation,
[`distrib_grad_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Weibull1Distrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
u <- (y / 2)^1.5
all.equal(g$mu, 1.5 * (u - 1) / 2)
#> [1] TRUE
all.equal(g$sigma, 1 / 1.5 + (1 - u) * log(y / 2))
#> [1] TRUE

# The summed score vanishes at the maximum likelihood estimate.
set.seed(5)
z <- distrib_rng(d, 500, th)
mle <- coef(fit_distrib(d, z))
vapply(distrib_gradient(d, z, as.list(mle)), sum, numeric(1))
#>            mu         sigma 
#> -2.907499e-07  3.387017e-06 

# On the link scale both components are multiplied by h' = theta, both
# parameters riding a log by default.
distrib_gradient(d, y, th, scale = "link")$mu / g$mu
#> [1] 2 2 2
```
