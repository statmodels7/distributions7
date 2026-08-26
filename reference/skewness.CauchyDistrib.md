# The Cauchy Distribution Has No Skewness

Returns `NaN`. Skewness is the third standardized central moment, so it
needs the first three moments of \\Y\\, and the Cauchy has none of them.
Its density decays like \\y^{-2}\\, so \\\int \|y\|^p
f(y)\\\mathrm{d}y\\ diverges for every \\p \ge 1\\; already the mean
fails to exist, and the second and third moments fail worse. `NaN` is
the value of a quantity that is not defined, and it is returned directly
so that no quadrature is attempted.

## Arguments

- x:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- theta:

  A named list of parameters with components `mu` and `sigma`, each a
  numeric vector of length 1 or `n`. The values are not read, only their
  lengths, but they are still aligned and validated, so a missing or
  out-of-bounds component throws exactly as it would for a family whose
  moments exist.

- ...:

  Unused, and accepted so that the signature matches the generic's.
  Arguments meaningful to
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
  have no effect here, no quadrature being run.

## Value

A numeric vector of `NaN`, of length
`max(length(theta$mu), length(theta$sigma))`, matching the
one-value-per- observation shape every moment method returns.

## Details

The generic
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
evaluates \\\gamma_1 = \mathbb{E}\[((Y -
\mathbb{E}Y)/\mathrm{sd}(Y))^3\]\\ for a `distrib` object by calling
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
three times. On a Cauchy each of those integrals diverges, and a
numerical quadrature over a divergent integral does not fail loudly: it
returns whatever its truncation happens to give, a number that changes
with the panel layout and looks like an estimate. This method
short-circuits that.

The same holds for
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md),
[`variance.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.CauchyDistrib.md)
and
[`kurtosis.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.CauchyDistrib.md),
and for the sample versions: a Cauchy sample mean does not converge as
\\n\\ grows, it is itself Cauchy with the same scale, so no amount of
data settles it. The example below shows that.

What the Cauchy does have is a **median**, equal to \\\mu\\, and a
half-interquartile range, equal to \\\sigma\\. Both are exact and finite
at every parameter value, and both are available from
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md).
Reach for those when a location and a spread are wanted.

## Notation

\\\mu\\ is the location and \\\sigma \> 0\\ the scale of the Cauchy, in
the parametrization
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
uses. Neither is a moment: \\\mu\\ is the median and \\\sigma\\ the
half-interquartile range.

## See also

[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md),
[`variance.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.CauchyDistrib.md)
and
[`kurtosis.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.CauchyDistrib.md),
which return `NaN` for the same reason;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the median and quartiles, which do exist;
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
for the generic and the families that answer it with a number.

## Examples

``` r
d <- cauchy_distrib()

# No moment of the Cauchy exists, so all four are NaN.
skewness(d, list(mu = 0, sigma = 1))
#> [1] NaN
c(mean(d, list(mu = 0, sigma = 1)),
  variance(d, list(mu = 0, sigma = 1)),
  kurtosis(d, list(mu = 0, sigma = 1)))
#> [1] NaN NaN NaN

# One value per observation, as for a family whose moments do exist.
skewness(d, list(mu = c(0, 1, 2), sigma = 1))
#> [1] NaN NaN NaN

# Why: a Cauchy sample mean is itself Cauchy, so it never settles.
set.seed(1)
y <- distrib_rng(d, 1e6, list(mu = 0, sigma = 1))
vapply(c(1e3, 1e4, 1e5, 1e6), function(n) mean(y[1:n]), numeric(1))
#> [1] -0.641804492  1.220622423 -2.147170390  0.006069028

# The median and the half-IQR are exact, and are mu and sigma.
distrib_quantile(d, 0.5, list(mu = 3, sigma = 2))
#> [1] 3
q <- distrib_quantile(d, c(0.25, 0.75), list(mu = 3, sigma = 2))
diff(q) / 2
#> [1] 2
```
