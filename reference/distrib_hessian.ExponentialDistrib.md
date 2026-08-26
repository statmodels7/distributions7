# Exponential Observed Hessian

Computes the second derivative of the exponential log-density with
respect to the mean, one value per observation, in closed form:
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\mu^2} -
\dfrac{2y}{\mu^3} = \dfrac{\mu - 2y}{\mu^3}.\$\$ It is positive wherever
\\y \< \mu/2\\, so a single observation below half the mean contributes
convexity; summed over a sample the curvature is negative at the
estimate, where \\\bar y = \mu\\. The expectation is
[`distrib_expected_hessian.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ExponentialDistrib.md).

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations, non-negative.

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

[`distrib_gradient.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ExponentialDistrib.md)
for the score,
[`distrib_expected_hessian.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ExponentialDistrib.md)
for the expectation of this quantity, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
y <- c(0.3, 1.1, 4.0)
th <- list(mu = 2)
h <- distrib_hessian(d, y, th)

# The closed form, written out.
all.equal(h$mu_mu, (2 - 2 * y) / 2^3)
#> [1] TRUE

# Positive below mu/2, negative above it.
data.frame(y = y, mu_mu = h$mu_mu, below_half = y < 2 / 2)
#>     y  mu_mu below_half
#> 1 0.3  0.175       TRUE
#> 2 1.1 -0.025      FALSE
#> 3 4.0 -0.750      FALSE

# A central difference of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 2 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 2 - eps))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#> [1] TRUE
```
