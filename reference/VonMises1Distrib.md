# von Mises Distribution Class

The S7 class of the von Mises family, the natural distribution for an
angle, parametrized by a mean direction \\\mu\\ and a concentration
\\\kappa \> 0\\, with density \\f(y) = e^{\kappa\cos(y-\mu)}/\\2\pi
I_0(\kappa)\\\\ on \\\[-\pi, \pi)\\. It inherits from
`continuous_distrib`; the nine methods listed below are registered on it
directly.

This is the first family in the package whose support has the topology
of a circle: the two ends of the interval are the same point, so the
density need not vanish at either.

Build one with
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
VonMises1Distrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0)
)
```

## Arguments

- distrib_name:

  A single character string specifying the name of the distribution
  (e.g., `"student t"`).

- dimension:

  A character string indicating the dimensionality (`"univariate"` or
  `"multivariate"`).

- bounds:

  A numeric vector of length 2 defining the overall support of the
  distribution `c(lower, upper)`.

- params:

  A character vector containing the names of the distribution parameters
  (e.g., `c("mu", "sigma")`).

- params_interpretation:

  A character vector (typically named) providing the statistical
  interpretation of each parameter (e.g., `c(mu = "location")`).

- n_params:

  A numeric value specifying the total number of parameters.

- params_bounds:

  A list of numeric vectors of length 2, specifying the valid
  mathematical domain for each individual parameter.

- link_params:

  A list of link function objects corresponding to each parameter,
  primarily used to map parameters to the unconstrained real line for
  optimization algorithms.

- params_smooth:

  An optional named logical vector flagging, for each parameter, whether
  the log-likelihood is differentiable with respect to it. Defaults to
  all `TRUE` (leave empty). Set an entry to `FALSE` for parameters at
  which the log-likelihood has a kink (e.g. the location of a Laplace
  distribution): the observed Hessian is then degenerate and the
  expected information must be obtained from the score variance rather
  than from \\-\mathbb{E}\[H\]\\ (see
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

## Value

An S7 object of class `VonMises1Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
they hold `"von mises1"`, `"univariate"`, `c(-pi, pi)`,
`c("mu", "kappa")`, the interpretations
`c(mu = "mean direction", kappa = "concentration")`, `2`, and the
domains \\(-\pi, \pi)\\ and \\(0, \infty)\\.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises1Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.VonMises1Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.VonMises1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises1Distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.VonMises1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises1Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.VonMises1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises1Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.VonMises1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.VonMises1Distrib.md).

The **quantile** comes from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md),
by root finding on the distribution function; the distribution function
itself is this class's own and is a Bessel series, not the parent's
quadrature. The four moments come from
[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md)
and its siblings, numerically, and are the ordinary moments of \\Y\\ as
a number rather than the circular ones.

## See also

[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
to build one;
[`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
for the same law parametrized by the mean resultant length;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the analogous family on the line;
[`distrib_cdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises1Distrib.md)
for the Bessel series.

## Examples

``` r
d <- vonmises1_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The support is the circle, written as a half-open interval.
d@bounds
#> [1] -3.141593  3.141593
d@params
#> [1] "mu"    "kappa"
d@params_interpretation
#>               mu            kappa 
#> "mean direction"  "concentration" 

# The direction rides a bounded link and the concentration a log.
vapply(d@link_params, function(l) l@link_name, character(1))
#>                                                     mu 
#> "bounded(lwr=-3.14159265358979, upr=3.14159265358979)" 
#>                                                  kappa 
#>                                                  "log" 
```
