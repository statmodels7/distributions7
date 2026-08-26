# Laplace Observed Hessian, Rate Parametrization

Computes the three distinct second derivatives of the Laplace
log-density with respect to \\\mu\\ and \\\lambda\\, one value per
observation, in closed form. Writing \\r = y - \mu\\,
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = 0, \qquad
\dfrac{\partial^2 \ell}{\partial \lambda^2} = -\dfrac{1}{\lambda^2},
\qquad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial \lambda} =
\mathrm{sign}(r).\$\$

The first entry is exactly zero for every observation, the log-density
being piecewise linear in \\\mu\\; see
[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md)
for what that costs. The curvature in \\\lambda\\ carries no data at
all, so in this parametrization the observed and expected values agree
in that entry.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `lambda` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `lambda_lambda` and
`mu_lambda`, each of length `length(y)`. `mu_mu` is a vector of zeros
and `lambda_lambda` is constant at \\-1/\lambda^2\\.

## See also

[`distrib_expected_hessian.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Laplace2Distrib.md)
for the information, whose location entry differs from this one;
[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md)
for the scale parametrization; and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, lambda = 2)
h <- distrib_hessian(d, y, th)

# Zero in the location, constant in the rate, a sign in the mixed entry.
h$mu_mu
#> [1] 0 0 0
h$lambda_lambda
#> [1] -0.25 -0.25 -0.25
all.equal(h$mu_lambda, sign(y - 0.4))
#> [1] TRUE

# The rate entry already equals its own expectation; the location one
# does not, and the expected page says why.
rbind(observed = vapply(h, function(v) v[1], numeric(1)),
      expected = vapply(distrib_expected_hessian(d, y, th),
                        function(v) v[1], numeric(1)))
#>          mu_mu lambda_lambda mu_lambda
#> observed     0         -0.25        -1
#> expected    -4         -0.25         0
```
