# Gamma Score in Mean and Dispersion

Computes the first derivatives of the gamma log-density with respect to
\\\mu\\ and \\\phi\\, one value per observation, in closed form. With
\\s = 1/\phi\\ the shape and \\z = y/\mu\\, \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^2}, \qquad \dfrac{\partial
\ell}{\partial \phi} = -s^2\left\\\log s + 1 - \psi(s) + \log z -
z\right\\,\$\$ with \\\psi\\ the digamma function. The first is the
score of a gamma generalized linear model, the residual divided by the
variance function \\\phi\mu^2\\.

The dispersion component is a difference of two quantities that agree to
leading order as \\s\\ grows: \\\log s - \psi(s)\\ and \\\log z - (z -
1)\\ each go to zero at the boundary this family tends towards. The
kernel computes each of them as a polygamma minus its own asymptote, so
the digits survive at large shape where the direct difference loses
them.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

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
\phi\mu^2\\. \\\psi\\ is the digamma function, \\\psi = (\log\Gamma)'\\.

## See also

[`distrib_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma1Distrib.md)
for their expectation,
[`distrib_grad_y.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gamma1Distrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)
g <- distrib_gradient(d, y, th)

# The generalized linear model score: residual over the variance function.
all.equal(g$mu, (y - 3) / (0.5 * 3^2))
#> [1] TRUE

# The dispersion component, written out with the digamma function.
s <- 1 / 0.5
z <- y / 3
all.equal(g$phi, -s^2 * (log(s) + 1 - digamma(s) + log(z) - z))
#> [1] TRUE

# The mean component vanishes at y = mu, whatever the dispersion.
distrib_gradient(d, 3, list(mu = 3, phi = c(0.1, 0.5, 2)))$mu
#> [1] 0 0 0

# Summed over a fitted sample the score is at the optimizer's tolerance.
set.seed(4)
zz <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, zz)
vapply(distrib_gradient(d, zz, as.list(coef(fit))), sum, numeric(1))
#>            mu           phi 
#> -7.048322e-14 -6.346140e-06 
```
