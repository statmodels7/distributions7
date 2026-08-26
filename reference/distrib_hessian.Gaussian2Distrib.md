# Gaussian Observed Hessian in Mean and Variance

Computes the three distinct second derivatives of the Gaussian
log-density with respect to \\\mu\\ and \\v = \sigma^2\\, one value per
observation, in closed form. With \\r = y - \mu\\, \$\$\ell^{(\mu\mu)} =
-\dfrac{1}{v}, \qquad \ell^{(\mu v)} = -\dfrac{r}{v^2}, \qquad
\ell^{(vv)} = \dfrac{1}{2v^2} - \dfrac{r^2}{v^3}.\$\$ Only the curvature
in \\\mu\\ is free of the data; the other two vary with the residual,
and their expectations are
[`distrib_expected_hessian.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian2Distrib.md).

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- y:

  A numeric vector of observations.

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

A named list of three numeric vectors, `mu_mu`, `mu_sigma2` and
`sigma2_sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`. The three name the
distinct entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density with respect
to parameters \\i\\ and \\j\\. Parenthesized superscripts name
derivatives; a subscript on \\\ell\\ never does.

## See also

[`distrib_gradient.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian2Distrib.md)
for the score,
[`distrib_expected_hessian.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian2Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian2Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gaussian2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, sigma2 = 4)
h <- distrib_hessian(d, y, th)

# The curvature in mu is constant at -1/v; the other two are not.
h$mu_mu
#> [1] -0.25 -0.25 -0.25
r <- y - 1
all.equal(h$mu_sigma2, -r / 4^2)
#> [1] TRUE
all.equal(h$sigma2_sigma2, 1 / (2 * 4^2) - r^2 / 4^3)
#> [1] TRUE

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-5
up <- distrib_gradient(d, y, list(mu = 1, sigma2 = 4 + eps))$sigma2
dn <- distrib_gradient(d, y, list(mu = 1, sigma2 = 4 - eps))$sigma2
all.equal((up - dn) / (2 * eps), h$sigma2_sigma2, tolerance = 1e-6)
#> [1] TRUE
```
