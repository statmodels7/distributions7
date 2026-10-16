# Pseudo-Huber Distribution Class, Standard-Deviation Parametrization

This is the S7 class of the pseudo-Huber family parametrized by its
location \\\mu\\, its standard deviation \\\sigma\\ and its shape
\\\nu\\. It is the family of
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
with the scale \\\sigma_1\\ of that family replaced by \\\sigma =
\sigma_1 \sqrt{R(\nu)}\\, where \\R(\nu) = \sqrt{\nu}\\ K_2(\sqrt{\nu})
/ K_1(\sqrt{\nu})\\. It inherits from `continuous_distrib`.

Build one with
[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md),
which supplies the three link functions and fills the properties in.
This page documents the raw S7 constructor, which takes the parent's
properties; the validator of `distrib` checks them, and the constructor
adds no check of its own.

## Usage

``` r
PseudoHuber2Distrib(
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

An S7 object of class `PseudoHuber2Distrib`, inheriting from
`continuous_distrib` and from `distrib`. For an object built by
[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md)
the properties hold `"pseudo huber2"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma", "nu")`, the interpretations
`c(mu = "location", sigma = "standard deviation", nu = "shape")`, `3`,
and the domains \\(-\infty, \infty)\\, \\(0, \infty)\\, \\(0, \infty)\\.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuber2Distrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuber2Distrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuber2Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PseudoHuber2Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuber2Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuber2Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuber2Distrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuber2Distrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuber2Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuber2Distrib.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuber2Distrib.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuber2Distrib.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuber2Distrib.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuber2Distrib.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuber2Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuber2Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuber2Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuber2Distrib.md);

- the fitting helpers:
  [`distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.PseudoHuberDistrib.md)
  and
  [`expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md)
to build one;
[PseudoHuberDistrib](https://statmodels7.github.io/distributions7/reference/PseudoHuberDistrib.md)
for the scale parametrization.

## Examples

``` r
d <- pseudohuber2_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE
d@params_interpretation
#>                   mu                sigma                   nu 
#>           "location" "standard deviation"              "shape" 
```
