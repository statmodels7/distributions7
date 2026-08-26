# Laplace Observed Hessian

Computes the three distinct second derivatives of the Laplace
log-density with respect to \\\mu\\ and \\\sigma\\, one value per
observation, in closed form. Writing \\r = y - \mu\\,
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = 0, \qquad
\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma -
2\|r\|}{\sigma^3}, \qquad \dfrac{\partial^2 \ell}{\partial \mu \\
\partial \sigma} = -\dfrac{\mathrm{sign}(r)}{\sigma^2}.\$\$

The first entry is **exactly zero for every observation**, and that is
the correct value: the log-density is piecewise linear in \\\mu\\, so
away from \\r = 0\\ its second derivative vanishes, and at \\r = 0\\ the
derivative does not exist. The curvature of the log-likelihood in the
location is concentrated in a set of measure zero, and no expectation of
this quantity sees it.

The consequence is that the observed information here is **not** the
information.
[`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md)
returns \\-1/\sigma^2\\, obtained from the variance of the score, and
that page explains which identity holds and which fails.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

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
`mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
`mu_mu` is a vector of zeros.

## See also

[`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md),
which returns the information and is not the expectation of this;
[`distrib_gradient.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LaplaceDistrib.md)
for the score;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
the estimation route this family needs; and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
h <- distrib_hessian(d, y, th)

# The curvature in the location is exactly zero at every observation.
h$mu_mu
#> [1] 0 0 0

# The other two, written out.
r <- y - 0.4
all.equal(h$sigma_sigma, (1.5 - 2 * abs(r)) / 1.5^3)
#> [1] TRUE
all.equal(h$mu_sigma, -sign(r) / 1.5^2)
#> [1] TRUE

# Averaging the observed curvature in mu gives 0, and the information is
# 1/sigma^2. The two disagree because this family is not regular in mu.
set.seed(12)
z <- distrib_rng(d, 1e5, th)
c(observed_mean = mean(distrib_hessian(d, z, th)$mu_mu),
  expected = distrib_expected_hessian(d, 0, th)$mu_mu)
#> observed_mean      expected 
#>     0.0000000    -0.4444444 
```
