# Lognormal Observed Hessian

Computes the three distinct second derivatives of the lognormal
log-density with respect to \\\mu\\ and \\\sigma^2\\, one value per
observation, in closed form. With \\r = \log y - \mu\\,
\$\$\ell^{(\mu\mu)} = -\dfrac{1}{\sigma^2}, \qquad \ell^{(\mu\sigma^2)}
= -\dfrac{r}{\sigma^4}, \qquad \ell^{(\sigma^2\sigma^2)} =
\dfrac{1}{2\sigma^4} - \dfrac{r^2}{\sigma^6}.\$\$ These are the
Gaussian's at \\\log y\\, the Jacobian of the log transformation
carrying no parameter. Only the curvature in \\\mu\\ is free of the
data; the other two vary with the residual on the log scale, and their
expectations are
[`distrib_expected_hessian.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Lognormal1Distrib.md).

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
`mu_sigma2`, in that order, each of length
`max(length(y), length(mu), length(sigma2))`. The three name the
distinct entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.
\\\sigma^2\\ is the variance of \\\log Y\\.

## See also

[`distrib_gradient.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal1Distrib.md)
for the score,
[`distrib_expected_hessian.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Lognormal1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Lognormal1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)
h <- distrib_hessian(d, y, th)

# The curvature in mu is constant at -1/sigma2; the other two are not.
h$mu_mu
#> [1] -2.777778 -2.777778 -2.777778
r <- log(y) - 0.5
all.equal(h$mu_sigma2, -r / 0.36^2)
#> [1] TRUE
all.equal(h$sigma2_sigma2, 1 / (2 * 0.36^2) - r^2 / 0.36^3)
#> [1] TRUE

# Identical to the Gaussian's at log y.
all.equal(h, distrib_hessian(gaussian2_distrib(), log(y), th))
#> [1] "Names: 2 string mismatches"                     
#> [2] "Component 2: Mean relative difference: 1.049141"
#> [3] "Component 3: Mean relative difference: 2.802071"

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 0.5, sigma2 = 0.36 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 0.5, sigma2 = 0.36 - eps))$mu
all.equal((up - dn) / (2 * eps), h$mu_sigma2, tolerance = 1e-5)
#> [1] TRUE
```
