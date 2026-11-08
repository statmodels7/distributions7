# Raw and Central Moments of a Distribution

Computes the raw moment \\E\[Y^p\]\\ or the central moment
\\E\[(Y-\mu)^p\]\\ of a distribution at one or more parameter settings.
The expectation is taken numerically by
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
which integrates over the support of a continuous family and sums the
series of a discrete one. The order \\p\\ may be fractional where \\(y -
c)^p\\ is real over the whole support, where \\c\\ is the centering
value (for example a raw moment of a family on the positive half-line);
the integrand is \\(y - c)^p\\ and not \\\|y - c\|^p\\. Every setting
supplied in `theta` is evaluated in one batched call, and the result
carries one value per setting.

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
  per parameter of `distrib`. Each component has length 1 or a common
  length \\n\\, and `p` and `mu` likewise, so that \\n\\ settings are
  evaluated in one call; other lengths signal an error. The list is
  aligned and validated by name, so a missing or out-of-bounds component
  signals an error.

- p:

  The order of the moment. A numeric vector of length 1 or of the number
  of settings, recycled against the components of `theta`. Defaults to
  1, the mean. A non-integer order is accepted where \\(y - c)^p\\ is
  real over the support. A negative order on a family whose support
  reaches zero gives a divergent integral, and the quadrature signals an
  error.

- central:

  Should the moment be taken about the mean? It is a single logical
  whose default `FALSE` gives the raw moment and needs no extra pass.
  `TRUE` costs one further evaluation unless `mu` is supplied.

- mu:

  The centering value used when `central = TRUE`. It is a numeric vector
  of length 1 or of the number of settings, or `NULL` (the default),
  which computes the mean numerically. Read only when `central = TRUE`.

- ...:

  Passed to
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
  and from there to the quadrature or the series. Names here must not
  collide with the names in `theta`.

## Value

A numeric vector of moments, of length equal to the longest of the
components of `theta` and of `p`. Where the integral or the series
diverges, or the integrand is not real, an error is signaled in place of
a value.

## What is computed

With `central = FALSE`, \$\$m_p = \int y^p f(y \mid \theta)\\
\mathrm{d}y \quad\text{or}\quad \sum_y y^p\\ f(y \mid \theta),\$\$ and
with `central = TRUE` the same integral or sum with \\y\\ replaced by
\\y - \mu\\. The centering value \\\mu\\ is the mean, obtained by a
first call at `p = 1`, unless it is supplied through the `mu` argument.
Supplying `mu = 0` with `central = TRUE` therefore returns the raw
moment, and supplying a fitted or a theoretical mean saves one pass.

## Accuracy and cost

A continuous family's moment is computed by quadrature and a discrete
family's by a series, so the result is accurate to the tolerance of
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
and is not exact. A closed form, where a family registers one, is
cheaper by about two orders of magnitude. Every univariate family except
the two von Mises registers closed forms for
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
and
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md),
so this function is reached through them only for those two.

## A moment that does not exist

When the integral or the series diverges, or the integrand is not real
over the support, this function signals an error: the quadrature does
not reach the requested accuracy (the mean of a Cauchy, a negative order
on a family whose support reaches zero, a fractional central moment of a
Gaussian) or the series does not converge (a fractional central moment
of a Poisson). An integral that only appears to converge can still
return a number that depends on the panel layout. A family whose moments
fail to exist therefore registers methods that return `NaN` directly, as
the Cauchy does through
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md),
and a family whose moments exist only above a threshold returns `Inf` or
`NaN` below it, as the Student t does.

## Notation

\\Y\\ is the response, \\f(y \mid \theta)\\ its density or mass
function, \\\theta\\ the parameter on its own scale, and \\\mu =
E\[Y\]\\ the mean.

## See also

[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
for the quadrature and the series on which this rests;
[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the standard moments, which use a family's closed form when one is
registered.

## Examples

``` r
d <- gaussian1_distrib()

# These calls compute the first four moments of a Gaussian, raw and central.
moment(d, list(mu = 2, sigma = 3), p = 1)                  # 2
#> [1] 2
moment(d, list(mu = 2, sigma = 3), p = 2)                  # mu^2 + sigma^2
#> [1] 13
moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE)  # sigma^2
#> [1] 9

# p is recycled against theta, so one call covers a grid of orders.
moment(d, list(mu = 0, sigma = 1), p = 1:4, central = TRUE)   # 0, 1, 0, 3
#> [1] -2.775558e-17  1.000000e+00  0.000000e+00  3.000000e+00

# Centering at zero recovers the raw moment.
all.equal(moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE, mu = 0),
          moment(d, list(mu = 2, sigma = 3), p = 2))
#> [1] TRUE

# The call returns one value per parameter setting.
moment(d, list(mu = c(0, 1, 2), sigma = 1), p = 2, central = TRUE)
#> [1] 1 1 1

# On a discrete family the expectation is an exact sum.
all.equal(moment(poisson_distrib(), list(mu = 3), p = 2, central = TRUE), 3)
#> [1] TRUE
```
