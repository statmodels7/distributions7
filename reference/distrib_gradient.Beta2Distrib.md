# Beta Score in the Shapes

Computes the first derivatives of the beta log-density with respect to
\\\alpha\\ and \\\beta\\, one value per observation, in closed form:
\$\$\dfrac{\partial \ell}{\partial \alpha} = \log y - \psi(\alpha) +
\psi(\alpha+\beta), \qquad \dfrac{\partial \ell}{\partial \beta} =
\log(1-y) - \psi(\beta) + \psi(\alpha+\beta),\$\$ with \\\psi\\ the
digamma function. The beta is an exponential family in the shapes with
sufficient statistics \\\log y\\ and \\\log(1-y)\\, so each component is
the corresponding statistic minus its expectation. That is also why
every derivative beyond this one is free of the response.

The value is computed in plain R, this family carrying no compiled
kernel, and \\\log(1-y)\\ is formed with
[`base::log1p()`](https://rdrr.io/r/base/Log.html) so that it stays
accurate at `y` near zero.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. An endpoint makes a
  logarithm infinite and the score non-finite.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `alpha` and `beta`, each of length
`max(length(y), length(alpha), length(beta))`.

## Notation

\\\ell\\ is the log-density of one observation and \\\alpha, \beta \>
0\\ the two shapes. \\\psi\\ is the digamma function, \\\psi =
(\log\Gamma)'\\.

## See also

[`distrib_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md)
for the second derivatives,
[`distrib_gradient.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta1Distrib.md)
for the same score in the mean and the precision,
[`beta2_higher()`](https://statmodels7.github.io/distributions7/reference/beta2_higher.md)
for the orders above, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out with the digamma function.
all.equal(g$alpha, log(y) - digamma(2) + digamma(7))
#> [1] TRUE
all.equal(g$beta, log1p(-y) - digamma(5) + digamma(7))
#> [1] TRUE

# Each component is a sufficient statistic minus its expectation, so the
# sample mean of log(y) matches psi(alpha) - psi(alpha + beta).
set.seed(8)
z <- distrib_rng(d, 2e5, th)
c(sample = mean(log(z)), theory = digamma(2) - digamma(7))
#>    sample    theory 
#> -1.449634 -1.450000 

# Summed over a fitted sample the score is at the optimizer's tolerance.
zz <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, zz)
vapply(distrib_gradient(d, zz, as.list(coef(fit))), sum, numeric(1))
#>        alpha         beta 
#> 1.903932e-06 2.565011e-07 
```
