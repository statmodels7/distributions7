# Beta Score in Mean and Precision

Computes the first derivatives of the beta log-density with respect to
\\\mu\\ and \\\phi\\, one value per observation, in closed form. With
\\\alpha = \mu\phi\\, \\\beta = (1-\mu)\phi\\ and \\\psi\\ the digamma
function, \$\$\dfrac{\partial \ell}{\partial \mu} = \phi\left\\
\log\dfrac{y}{1-y} - \psi(\alpha) + \psi(\beta)\right\\,\$\$
\$\$\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) -
\mu\psi(\alpha) - (1-\mu)\psi(\beta) + \mu\log y + (1-\mu)\log(1-y).\$\$
The data enter only through \\\log y\\ and \\\log(1-y)\\, the beta being
an exponential family in the shapes, and the mean component sees them
only through the log-odds \\\log\\y/(1-y)\\\\.

The two expectations that make the score have mean zero are
\\\mathbb{E}\[\log Y\] = \psi(\alpha) - \psi(\phi)\\ and
\\\mathbb{E}\[\log(1-Y)\] = \psi(\beta) - \psi(\phi)\\.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. An endpoint makes a
  logarithm infinite and the score non-finite.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

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

\\\ell\\ is the log-density of one observation, \\\mu \in (0,1)\\ the
mean and \\\phi \> 0\\ the precision. \\\psi\\ is the digamma function,
\\\psi = (\log\Gamma)'\\.

## See also

[`distrib_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta1Distrib.md)
for their expectation,
[`distrib_grad_y.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Beta1Distrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out with the digamma function.
a <- 0.4 * 5
b <- 0.6 * 5
all.equal(g$mu, 5 * (log(y / (1 - y)) - digamma(a) + digamma(b)))
#> [1] TRUE
all.equal(g$phi, digamma(5) - 0.4 * digamma(a) - 0.6 * digamma(b) +
                 0.4 * log(y) + 0.6 * log(1 - y))
#> [1] TRUE

# The mean component vanishes where the log-odds equal psi(a) - psi(b).
distrib_gradient(d, plogis(digamma(a) - digamma(b)), th)$mu
#> [1] -1.110223e-15

# Summed over a fitted sample the score is at the optimizer's tolerance.
set.seed(6)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
#>           mu          phi 
#> 2.217664e-07 3.741311e-08 
```
