# NB1 Observed Hessian

Computes the three distinct second derivatives of the NB1 log-mass with
respect to \\\mu\\ and \\\theta\\, one value per observation, in closed
form. They are the two-variable chain rule through the size \\r =
\mu/\theta\\ applied once more: the second derivative of the log-mass in
\\r\\ is \\\psi'(y+r) - \psi'(r)\\, with \\\psi'\\ the trigamma
function, and the term in which \\\theta\\ appears outside \\r\\
contributes \\-1/(1+\theta)\\ to the mixed entry.

Because \\\partial r/\partial\mu = 1/\theta\\ and \\\partial
r/\partial\theta = -r/\theta\\, every component divides by a power of
\\\theta\\, and the components in \\\theta\\ lose their digits at small
\\\theta\\ for the reason given on
[`distrib_gradient.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin1Distrib.md).

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- y:

  A numeric vector of counts.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
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

A named list of three numeric vectors, `mu_mu`, `mu_theta` and
`theta_theta`, in that order, each of length
`max(length(y), length(mu), length(theta))`. The three name the distinct
entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-mass in parameters
\\i\\ and \\j\\; parenthesized superscripts name derivatives. \\r =
\mu/\theta\\ is the size, \\\psi\\ the digamma function and \\\psi'\\
the trigamma.

## See also

[`distrib_gradient.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin1Distrib.md)
for the score,
[`distrib_expected_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 4)
h <- distrib_hessian(d, y, th)
h
#> $mu_mu
#> [1]  0.00000000 -0.07812500 -0.09321181
#> 
#> $mu_theta
#> [1]  0.050589870  0.034964870 -0.009323325
#> 
#> $theta_theta
#> [1] -0.061179739  0.003195261  0.016858455
#> 

# The curvature in mu is exactly zero at y = 0, the size being 1 there, so
# the observed information is singular at that count.
h$mu_mu
#> [1]  0.00000000 -0.07812500 -0.09321181

# It is the second derivative of the log-mass, so a central difference of
# the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 4 + eps, theta = 4))$mu
dn <- distrib_gradient(d, y, list(mu = 4 - eps, theta = 4))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
#> [1] TRUE

# And the mixed entry against a difference of the other component.
up <- distrib_gradient(d, y, list(mu = 4, theta = 4 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 4, theta = 4 - eps))$mu
all.equal((up - dn) / (2 * eps), h$mu_theta, tolerance = 1e-5)
#> [1] TRUE
```
