# Gumbel Observed Hessian

Computes the three distinct second derivatives of the Gumbel log-density
with respect to \\\mu\\ and \\\sigma\\, one value per observation, in
closed form. With \\z = (y-\mu)/\sigma\\ and \\w = e^{-z}\\,
\$\$\ell^{(\mu\mu)} = -\dfrac{w}{\sigma^2}, \qquad \ell^{(\mu\sigma)} =
-\dfrac{1 - w + zw}{\sigma^2}, \qquad \ell^{(\sigma\sigma)} = \dfrac{1 -
2z + 2zw - z^2 w}{\sigma^2}.\$\$ The curvature in the location is
negative at every observation, \\w\\ being positive, so the
log-likelihood is concave in \\\mu\\ at any data set. The curvature in
the scale is not: it is positive wherever \\1 - 2z + 2zw - z^2w \> 0\\,
which includes \\y = \mu\\.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, in that order, each of length
`max(length(y), length(mu), length(sigma))`. The three name the distinct
entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.
\\z = (y-\mu)/\sigma\\ and \\w = e^{-z}\\.

## See also

[`distrib_gradient.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GumbelDistrib.md)
for the score,
[`distrib_expected_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GumbelDistrib.md)
for the expectation of this quantity,
[`distrib_deriv3.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GumbelDistrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)
h <- distrib_hessian(d, y, th)
h
#> $mu_mu
#> [1] -2.7182818 -1.0000000 -0.3678794
#> 
#> $sigma_sigma
#> [1] -5.1548455  1.0000000 -0.6321206
#> 
#> $mu_sigma
#> [1]  4.436564  0.000000 -1.000000
#> 

# The three closed forms, written out.
z <- (y - 0) / 1
w <- exp(-z)
all.equal(h$mu_mu, -w)
#> [1] TRUE
all.equal(h$mu_sigma, -(1 - w + z * w))
#> [1] TRUE
all.equal(h$sigma_sigma, 1 - 2 * z + 2 * z * w - z^2 * w)
#> [1] TRUE

# Concave in mu everywhere, and positive in sigma at y = mu.
c(all_negative_in_mu = all(h$mu_mu < 0),
  sigma_at_mu = distrib_hessian(d, 0, th)$sigma_sigma)
#> all_negative_in_mu        sigma_at_mu 
#>                  1                  1 

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 0 + eps, sigma = 1))$mu
dn <- distrib_gradient(d, y, list(mu = 0 - eps, sigma = 1))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
