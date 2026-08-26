# Beta-Binomial Observed Hessian in Its Shapes

Computes the three distinct second derivatives of the beta-binomial
log-mass with respect to the two shapes, one value per observation, in
closed form. They are the score's digamma differences with the trigamma
function in its place: \$\$\dfrac{\partial^2\ell}{\partial\alpha^2} =
\psi_1(y+\alpha) - \psi_1(\alpha) - \psi_1(n+S) + \psi_1(S), \qquad S =
\alpha+\beta,\$\$ and the same with \\n-y\\ and \\\beta\\. The **mixed**
component keeps only the shared term, \\-\psi_1(n+S) + \psi_1(S)\\, the
two shapes entering the log-mass separately otherwise, so it does not
depend on the data at all and equals its own expectation at every
observation.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

- y:

  A numeric vector of counts in \\\\0, \dots, n\\\\.

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

A named list of three numeric vectors, `alpha_alpha`, `beta_beta` and
`alpha_beta`, each of length
`max(length(y), length(alpha), length(beta))`. The three name the
distinct entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell\\ is the log-mass of one observation, \\\alpha, \beta \> 0\\ the
two beta shapes, \\n\\ the trial count and \\\psi_1\\ the trigamma
function, [`trigamma()`](https://rdrr.io/r/base/Special.html) in R.

## See also

[`distrib_gradient.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom2Distrib.md)
for the score,
[`distrib_expected_hessian.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaBinom2Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom2Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
th <- list(alpha = 2, beta = 3)
h <- distrib_hessian(d, 0:10, th)
names(h)
#> [1] "alpha_alpha" "beta_beta"   "alpha_beta" 

# The mixed component is free of the data: one number, repeated.
unique(h$alpha_beta)
#> [1] 0.1523847
-trigamma(10 + 5) + trigamma(5)
#> [1] 0.1523847

# A central difference of the score reproduces the pure-alpha component.
eps <- 1e-5
up <- distrib_gradient(d, 0:10, list(alpha = 2 + eps, beta = 3))$alpha
dn <- distrib_gradient(d, 0:10, list(alpha = 2 - eps, beta = 3))$alpha
all.equal((up - dn) / (2 * eps), h$alpha_alpha, tolerance = 1e-6)
#> [1] TRUE
```
