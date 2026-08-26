# Gamma Score in Mean and Variance

Computes the first derivatives of the gamma log-density with respect to
\\\mu\\ and \\\sigma^2\\, one value per observation, in closed form.
With \\\alpha = \mu^2/\sigma^2\\ and \\\lambda = \mu/\sigma^2\\,
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{-2\mu\psi(\alpha) +
2\mu\log\lambda + \mu + 2\mu\log y - y} {\sigma^2},\$\$
\$\$\dfrac{\partial \ell}{\partial \sigma^2} =
-\dfrac{\mu\left\\-\mu\psi(\alpha) + \mu + \mu(\log\lambda + \log y) -
y\right\\}{(\sigma^2)^2},\$\$ with \\\psi\\ the digamma function. Both
components carry the digamma, because both parameters move the shape; in
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md),
where the second parameter is a dispersion, only the second component
does.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
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

A named list of two numeric vectors, `mu` and `sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the mean
and \\\sigma^2 \> 0\\ the variance. \\\psi\\ is the digamma function,
\\\psi = (\log\Gamma)'\\.

## See also

[`distrib_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma2Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma2Distrib.md)
for their expectation,
[`distrib_gradient.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma1Distrib.md)
for the same score in the dispersion, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out with the digamma function.
al <- 9 / 2
lam <- 3 / 2
all.equal(g$mu,
          (-2 * 3 * digamma(al) + 2 * 3 * log(lam) + 3 +
             2 * 3 * log(y) - y) / 2)
#> [1] TRUE
all.equal(g$sigma2,
          -3 * (-3 * digamma(al) + 3 + 3 * (log(lam) + log(y)) - y) / 2^2)
#> [1] TRUE

# Summed over a fitted sample the score is at the optimizer's tolerance.
set.seed(5)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
#>            mu        sigma2 
#>  1.200567e-09 -8.849422e-10 
```
