# Raw and Central Moments of a Distribution

Computes the raw moment \\E\[Y^p\]\\ or the central moment
\\E\[(Y-\mu)^p\]\\ of a distribution at one or more parameter settings.
The expectation is taken numerically by
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
which integrates over the support of a continuous family and sums the
series of a discrete one. The order \\p\\ need not be a whole number, so
fractional and absolute moments are reachable from the same function.
Every setting supplied in `theta` is evaluated in one batched call, and
the result carries one value per setting.

## Usage

``` r
moment(distrib, theta, p = 1, central = FALSE, mu = NULL, ...)
```

## Arguments

- distrib:

  An object inheriting from `distrib`, from any of the family
  constructors such as
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
  or
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- theta:

  A named list of parameters on the parameter scale, with one component
  per parameter of `distrib`. Components may be vectors and are recycled
  against one another and against `p`, so several settings are evaluated
  in one call. The list is aligned and validated by name, so a missing
  or out-of-bounds component throws.

- p:

  The order of the moment. A numeric vector of length 1 or of the number
  of settings, recycled against the components of `theta`. Defaults to
  1, the mean. Non-integer orders are accepted; a negative order is
  accepted too and diverges for any family whose support reaches zero.

- central:

  Should the moment be taken about the mean? A single logical, `FALSE`
  by default, which gives the raw moment and needs no extra pass. `TRUE`
  costs one further evaluation unless `mu` is supplied.

- mu:

  The centering value used when `central = TRUE`. A numeric vector of
  length 1 or of the number of settings, or `NULL` (the default), which
  computes the mean numerically. Read only when `central = TRUE`.

- ...:

  Passed to
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
  and from there to the quadrature or the series. Names here must not
  collide with the names in `theta`.

## Value

A numeric vector of moments, of length equal to the longest of the
components of `theta` and of `p`. `NaN` where the integrand is not
defined, `Inf` where the integral diverges to infinity.

## What is computed

With `central = FALSE`, \$\$m_p = \int y^p f(y \mid \theta)\\
\mathrm{d}y \quad\text{or}\quad \sum_y y^p\\ f(y \mid \theta),\$\$ and
with `central = TRUE` the same integral or sum with \\y\\ replaced by
\\y - \mu\\. The centering value \\\mu\\ is the mean, obtained by a
first call at `p = 1`, unless it is supplied through the `mu` argument.
Supplying `mu = 0` with `central = TRUE` therefore returns the raw
moment, and supplying a fitted or a theoretical mean saves one pass.

## Accuracy and cost

A continuous family's moment is a quadrature and a discrete family's is
a series, so neither is exact. On a Gaussian the second and fourth
central moments come back at \\9 - 2.8\times10^{-14}\\ and \\3 -
8.9\times10^{-15}\\ against the exact 9 and 3; on a Poisson the series
is exact to the last bit at ordinary means. The price is the evaluation:
one numerical variance of a Gaussian costs 2.3 ms against 24
microseconds for the closed form the family registers, a factor of about
94, and one numerical skewness costs 3.5 ms against the same 24
microseconds. That gap is why 43 of the 45 shipped families answer
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
and
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
with a formula of their own; the two von Mises families are the ones
that reach this function.

## A moment that does not exist

A divergent integral does not announce itself. The quadrature returns
whatever its truncation gives, a number that moves with the panel layout
and looks like an estimate. A family whose moments fail to exist
therefore registers a method that returns `NaN` directly, as the Cauchy
does through
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md),
and a family whose moments exist only above a threshold returns `Inf`
where they do not, as the Student t does at \\\nu \le 4\\ for the
kurtosis. Reading a number back from this function on a family with no
analytical method is a statement about the quadrature and not about the
law.

## Notation

\\Y\\ is the response, \\f(y \mid \theta)\\ its density or mass
function, \\\theta\\ the parameter on its own scale, and \\\mu =
E\[Y\]\\ the mean.

## See also

[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
for the quadrature and the series this rests on;
[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the four standard moments, each of which prefers a family's closed
form when one is registered.

## Examples

``` r
d <- gaussian1_distrib()

# The first four moments of a Gaussian, raw and central.
moment(d, list(mu = 2, sigma = 3), p = 1)                  # 2
#> [1] 2
moment(d, list(mu = 2, sigma = 3), p = 2)                  # mu^2 + sigma^2
#> [1] 13
moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE)  # sigma^2
#> [1] 9

# p is recycled against theta, so one call covers a grid of orders.
moment(d, list(mu = 0, sigma = 1), p = 1:4, central = TRUE)   # 0, 1, 0, 3
#> [1] 0 1 0 3

# Centering at zero recovers the raw moment.
all.equal(moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE, mu = 0),
          moment(d, list(mu = 2, sigma = 3), p = 2))
#> [1] TRUE

# One value per parameter setting.
moment(d, list(mu = c(0, 1, 2), sigma = 1), p = 2, central = TRUE)
#> [1] 1 1 1

# On a discrete family the expectation is an exact sum.
all.equal(moment(poisson_distrib(), list(mu = 3), p = 2, central = TRUE), 3)
#> [1] TRUE
```
