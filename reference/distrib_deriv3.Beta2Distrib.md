# Beta Third-Order Derivatives in the Shapes

Computes the four distinct third derivatives of the beta log-density
with respect to \\\alpha\\ and \\\beta\\, in closed form:
\$\$\ell^{(\alpha\alpha\alpha)} = \psi''(\alpha+\beta) - \psi''(\alpha),
\qquad \ell^{(\alpha\alpha\beta)} = \ell^{(\alpha\beta\beta)} =
\psi''(\alpha+\beta), \qquad \ell^{(\beta\beta\beta)} =
\psi''(\alpha+\beta) - \psi''(\beta),\$\$ with \\\psi''\\ the second
derivative of the digamma function. The two mixed components are equal,
every mixed derivative of \\\log B(\alpha,\beta)\\ of a given order
being the same polygamma of the sum.

The values are free of the response, so `expected` selects nothing: the
same computation runs either way and the two results are identical to
the bit. `approx` and `nsim` are ignored for the same reason.

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
  expected third derivatives being the same numbers. Defaults to
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

A named list of four numeric vectors, `alpha_alpha_alpha`,
`alpha_alpha_beta`, `alpha_beta_beta` and `beta_beta_beta`, each of
length `length(y)` and constant within itself when the parameters are.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. \\\psi\\ is the digamma function and
\\\psi^{(m)}\\ its \\m\\th derivative.

## See also

[`distrib_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md)
for the order below,
[`distrib_deriv4.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Beta2Distrib.md)
for the order above,
[`beta2_higher()`](https://statmodels7.github.io/distributions7/reference/beta2_higher.md),
which computes this, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)
d3 <- distrib_deriv3(d, y, th)

# Four constants, and the two mixed components are the same number.
lapply(d3, unique)
#> $alpha_alpha_alpha
#> [1] 0.3805833
#> 
#> $alpha_alpha_beta
#> [1] -0.02353047
#> 
#> $alpha_beta_beta
#> [1] -0.02353047
#> 
#> $beta_beta_beta
#> [1] 0.02525926
#> 
psigamma(7, deriv = 2)
#> [1] -0.02353047

# Free of the response, so asking for the expectation changes nothing.
identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(alpha = 2 + eps, beta = 5))$alpha_alpha
dn <- distrib_hessian(d, y, list(alpha = 2 - eps, beta = 5))$alpha_alpha
all.equal((up - dn) / (2 * eps), d3$alpha_alpha_alpha, tolerance = 1e-4)
#> [1] TRUE
```
