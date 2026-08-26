# Geometric Observed Hessian

Computes the second derivative of the geometric log-mass with respect to
the mean, one value per observation, in closed form. Differentiating
\\(y-\mu)/(\mu(1+\mu))\\ gives \$\$\dfrac{\partial^2 \ell}{\partial
\mu^2} = \dfrac{-\mu(1+\mu) - (y-\mu)(1+2\mu)}{\mu^2(1+\mu)^2},\$\$
which depends on the data through \\y\\ alone. Its expectation is
[`distrib_expected_hessian.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GeometricDistrib.md).

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- y:

  A numeric vector of counts.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu`, of length
`max(length(y), length(mu))`.

## See also

[`distrib_gradient.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GeometricDistrib.md)
for the score,
[`distrib_expected_hessian.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GeometricDistrib.md)
for the expectation of this quantity, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- geometric_distrib()
y <- c(0, 2, 7)
th <- list(mu = 3)

# A central difference of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 3 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 3 - eps))$mu
all.equal((up - dn) / (2 * eps), distrib_hessian(d, y, th)$mu_mu,
          tolerance = 1e-6)
#> [1] TRUE

# It varies with the count, unlike the expected value.
rbind(observed = distrib_hessian(d, y, th)$mu_mu,
      expected = distrib_expected_hessian(d, y, th)$mu_mu)
#>                 [,1]        [,2]        [,3]
#> observed  0.06250000 -0.03472222 -0.27777778
#> expected -0.08333333 -0.08333333 -0.08333333
```
