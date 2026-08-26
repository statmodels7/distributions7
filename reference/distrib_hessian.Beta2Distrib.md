# Beta Observed Hessian in the Shapes

Computes the three distinct second derivatives of the beta log-density
with respect to \\\alpha\\ and \\\beta\\, in closed form:
\$\$\ell^{(\alpha\alpha)} = \psi'(\alpha+\beta) - \psi'(\alpha), \qquad
\ell^{(\alpha\beta)} = \psi'(\alpha+\beta), \qquad \ell^{(\beta\beta)} =
\psi'(\alpha+\beta) - \psi'(\beta),\$\$ with \\\psi'\\ the trigamma
function. **All three are free of the response**, the data entering the
log-density only through a term linear in the parameters. The values are
constant within a parameter setting and are recycled to the length of
`y`, and they equal
[`distrib_expected_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta2Distrib.md)
exactly.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `alpha_alpha`, `beta_beta` and
`alpha_beta`, in that order, each of length `length(y)`. The three name
the distinct entries of a symmetric \\2 \times 2\\ matrix per
observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.
\\\psi'\\ is the trigamma function.

## See also

[`distrib_gradient.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta2Distrib.md)
for the score, which is the one quantity here that does read the data;
[`distrib_expected_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta2Distrib.md),
which returns the same numbers;
[`beta2_higher()`](https://statmodels7.github.io/distributions7/reference/beta2_higher.md),
which computes this; and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)
h <- distrib_hessian(d, y, th)
h
#> $alpha_alpha
#> [1] -0.4913889 -0.4913889 -0.4913889
#> 
#> $beta_beta
#> [1] -0.06777778 -0.06777778 -0.06777778
#> 
#> $alpha_beta
#> [1] 0.1535452 0.1535452 0.1535452
#> 

# Written out with the trigamma function.
c(trigamma(7) - trigamma(2), trigamma(7) - trigamma(5), trigamma(7))
#> [1] -0.49138889 -0.06777778  0.15354518

# Free of the response, so it equals its own expectation to the bit.
identical(h, distrib_expected_hessian(d, y, th))
#> [1] TRUE

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(alpha = 2 + eps, beta = 5))$alpha
dn <- distrib_gradient(d, y, list(alpha = 2 - eps, beta = 5))$alpha
all.equal((up - dn) / (2 * eps), h$alpha_alpha, tolerance = 1e-5)
#> [1] TRUE
```
