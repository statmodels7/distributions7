# Beta Fourth-Order Derivatives in the Shapes

Computes the five distinct fourth derivatives of the beta log-density
with respect to \\\alpha\\ and \\\beta\\, in closed form:
\$\$\ell^{(\alpha^4)} = \psi'''(\alpha+\beta) - \psi'''(\alpha), \qquad
\ell^{(\beta^4)} = \psi'''(\alpha+\beta) - \psi'''(\beta),\$\$ and every
one of the three mixed components equal to \\\psi'''(\alpha+\beta)\\,
with \\\psi'''\\ the third derivative of the digamma function.

As at third order the values are free of the response, so `expected`,
`approx` and `nsim` are all without effect and the result is identical
to the bit whichever is asked for.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- expected:

  Logical of length 1, and without effect here, the observed and
  expected fourth derivatives being the same numbers. Defaults to
  `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of five numeric vectors, `alpha_alpha_alpha_alpha`,
`alpha_alpha_alpha_beta`, `alpha_alpha_beta_beta`,
`alpha_beta_beta_beta` and `beta_beta_beta_beta`, each of length
`length(y)` and constant within itself when the parameters are.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. \\\psi\\ is the digamma
function and \\\psi^{(m)}\\ its \\m\\th derivative.

## See also

[`distrib_deriv3.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Beta2Distrib.md)
for the order below,
[`beta2_higher()`](https://statmodels7.github.io/distributions7/reference/beta2_higher.md),
which computes this, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)
d4 <- distrib_deriv4(d, y, th)

# Five constants, and the three mixed components are one number.
lapply(d4, unique)
#> $alpha_alpha_alpha_alpha
#> [1] -0.4867412
#> 
#> $alpha_alpha_alpha_beta
#> [1] 0.007198199
#> 
#> $alpha_alpha_beta_beta
#> [1] 0.007198199
#> 
#> $alpha_beta_beta_beta
#> [1] 0.007198199
#> 
#> $beta_beta_beta_beta
#> [1] -0.01422963
#> 
psigamma(7, deriv = 3)
#> [1] 0.007198199

# Free of the response, so asking for the expectation changes nothing.
identical(d4, distrib_deriv4(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(alpha = 2 + eps, beta = 5))$alpha_alpha_alpha
dn <- distrib_deriv3(d, y, list(alpha = 2 - eps, beta = 5))$alpha_alpha_alpha
all.equal((up - dn) / (2 * eps), d4$alpha_alpha_alpha_alpha,
          tolerance = 1e-3)
#> [1] TRUE
```
