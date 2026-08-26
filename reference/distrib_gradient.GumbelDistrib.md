# Gumbel Score

Computes the first derivatives of the Gumbel log-density with respect to
\\\mu\\ and \\\sigma\\, one value per observation, in closed form. With
\\z = (y-\mu)/\sigma\\ and \\w = e^{-z}\\, \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{1 - w}{\sigma}, \qquad \dfrac{\partial
\ell}{\partial \sigma} = \dfrac{z(1 - w) - 1}{\sigma}.\$\$

The whole family is written in \\w\\, and under the model \\w\\ is
**standard exponential** whatever the parameters. Every expectation this
family needs is therefore a derivative of \\\Gamma\\ at 2, and the same
fact bounds the score above in \\\mu\\ and leaves it unbounded below: a
value far to the left makes \\w\\ enormous.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

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

## Value

A named list of two numeric vectors, `mu` and `sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\sigma \> 0\\ the scale. \\z = (y-\mu)/\sigma\\ is the standardized
value and \\w = e^{-z}\\, which is standard exponential under the model.

## See also

[`distrib_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GumbelDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GumbelDistrib.md)
for their expectation,
[`distrib_grad_y.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GumbelDistrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
z <- (y - 0) / 1
w <- exp(-z)
all.equal(g$mu, (1 - w) / 1)
#> [1] TRUE
all.equal(g$sigma, (z * (1 - w) - 1) / 1)
#> [1] TRUE

# The location component vanishes at y = mu, where w = 1.
distrib_gradient(d, 0, th)$mu
#> [1] 0

# It is bounded above by 1/sigma and unbounded below.
distrib_gradient(d, c(-3, 10), th)$mu
#> [1] -19.0855369   0.9999546

# Summed over a fitted sample the score is at the optimizer's tolerance.
set.seed(8)
s <- distrib_rng(d, 2000, list(mu = 3, sigma = 2))
fit <- fit_distrib(d, s)
vapply(distrib_gradient(d, s, as.list(coef(fit))), sum, numeric(1))
#>            mu         sigma 
#>  5.607183e-06 -2.795275e-05 
```
