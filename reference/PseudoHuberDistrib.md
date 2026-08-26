# Pseudo-Huber Distribution Class

The S7 class of the pseudo-Huber family, a location-scale family on the
whole real line whose log-density is the negative pseudo-Huber loss
\\-\sqrt{\nu + z^2}\\, with a shape \\\nu \> 0\\ interpolating between a
Laplace at \\\nu \to 0\\ and a Gaussian at \\\nu \to \infty\\. It is the
symmetric hyperbolic distribution, a special case of the generalized
hyperbolic. It inherits from `continuous_distrib`; the eleven methods
listed below are registered on it in this file.

Build one with
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md),
which supplies the three link functions and fills the properties in.
This page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
PseudoHuberDistrib(
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

An S7 object of class `PseudoHuberDistrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
they hold `"pseudo huber"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma", "nu")`, the interpretations
`c(mu = "location", sigma = "scale", nu = "shape")`, `3`, and the
domains \\(-\infty, \infty)\\, \\(0, \infty)\\, \\(0, \infty)\\.

## Methods

Registered on this class in this file:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuberDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PseudoHuberDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.PseudoHuberDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PseudoHuberDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuberDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuberDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.PseudoHuberDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuberDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuberDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PseudoHuberDistrib.md),
and the predicate
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.PseudoHuberDistrib.md),
which answers `FALSE` here.

Registered from other files: the mixed derivative
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.PseudoHuberDistrib.md)
and the distribution-function derivatives
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuberDistrib.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.PseudoHuberDistrib.md),
plus the four moments
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuberDistrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.PseudoHuberDistrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.PseudoHuberDistrib.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PseudoHuberDistrib.md)
in `moments.R`.

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
to build one;
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the two limits;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
for the other robust three-parameter family here;
[`distrib_pdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuberDistrib.md)
for the density.

## Examples

``` r
d <- pseudohuber_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"    "sigma" "nu"   
d@params_interpretation
#>         mu      sigma         nu 
#> "location"    "scale"    "shape" 

# This is the one family here whose expected information is not written
# out, and the predicate says so.
distributions7:::expected_hessian_exact(d)
#> [1] FALSE
```
