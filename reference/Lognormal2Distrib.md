# Lognormal Distribution Class, Mean and Variance of Y

The S7 class of the lognormal family parametrized by the mean \\m\\ and
the variance \\v\\ of \\Y\\. It inherits from `continuous_distrib`.
Build one with
[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md);
this page documents the raw S7 constructor; the validator of `distrib`
checks that the properties are consistent (the parameter count, a bound
and a link for every parameter, ordered support bounds).

## Usage

``` r
Lognormal2Distrib(
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

An S7 object of class `Lognormal2Distrib`, inheriting from
`continuous_distrib`. For an object built by
[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
the properties hold `"lognormal2"`, `"univariate"`, `c(0, Inf)`,
`c("mean", "var")`, the interpretations
`c(mean = "mean", var = "variance")`, `2` and the domains \\(0,
\infty)\\ for both.

## Methods

Registered on this class, each linked to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal2Distrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal2Distrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal2Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal2Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal2Distrib.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Lognormal2Distrib.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Lognormal2Distrib.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Lognormal2Distrib.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Lognormal2Distrib.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.Lognormal2Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/mean.Lognormal2Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/mean.Lognormal2Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/mean.Lognormal2Distrib.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
to build one;
[Lognormal1Distrib](https://statmodels7.github.io/distributions7/reference/Lognormal1Distrib.md)
for the parametrization on the log scale.

## Examples

``` r
d <- lognormal2_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE
d@params
#> [1] "mean" "var" 
```
