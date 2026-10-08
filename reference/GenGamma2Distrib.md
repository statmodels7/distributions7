# Generalized Gamma Distribution Class, Mean

The S7 class of the generalized gamma family parametrized by its mean
\\m\\, its shape \\d\\ and its power \\p\\. It inherits from
`continuous_distrib`. Build one with
[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md);
this page documents the raw S7 constructor, which validates none of the
relationships between the properties.

## Usage

``` r
GenGamma2Distrib(
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

An S7 object of class `GenGamma2Distrib`, inheriting from
`continuous_distrib`. For an object built by
[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md)
the properties hold `"gengamma2"`, `"univariate"`, `c(0, Inf)`,
`c("mean", "d", "p")`, the interpretations
`c(mean = "mean", d = "shape", p = "power")`, `3` and the domain \\(0,
\infty)\\ for each parameter.

## Methods

Registered on this class in this file:
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
the derivatives of the expected information through
[`register_dexpected()`](https://statmodels7.github.io/distributions7/reference/register_dexpected.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md),
[`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md),
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
and the four moments [`mean()`](https://rdrr.io/r/base/mean.html),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md).
The derivatives of the distribution function are the numerical ones of
the base class, as for
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

## See also

[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md)
to build one;
[GenGamma1Distrib](https://statmodels7.github.io/distributions7/reference/GenGamma1Distrib.md)
for the parametrization by the scale.

## Examples

``` r
d <- gengamma2_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE
d@params
#> [1] "mean" "d"    "p"   
```
