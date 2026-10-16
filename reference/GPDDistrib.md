# Generalized Pareto Distribution Class

This is the S7 class of the generalized Pareto family, the law of
exceedances over a high threshold, parametrized by a scale \\\sigma \>
0\\ and a shape \\\xi\\: \\f(y) = \sigma^{-1}(1 + \xi
y/\sigma)^{-1/\xi-1}\\ on \\y \ge 0\\. At \\\xi = 0\\ it is the
exponential, reached by a series, so a fit may move the shape through
zero.

For \\\xi \< 0\\ the support ends at \\-\sigma/\xi\\, so it moves with
the parameters. Because of this the family is non-regular, and for the
same reason its expected information exists only above \\\xi = -1/2\\.

Build one with
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md),
which supplies the two link functions. This page documents the raw S7
constructor, which takes the parent's properties; the validator of
`distrib` checks them, and the constructor adds no check of its own.

## Usage

``` r
GPDDistrib(
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

An S7 object of class `GPDDistrib`, inheriting from `continuous_distrib`
and from `distrib`. For an object built by
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
the properties hold `"generalized pareto"`, `"univariate"`, `c(0, Inf)`,
`c("sigma", "xi")`, the interpretations
`c(sigma = "scale", xi = "shape")`, `2`, and the domains \\(0,\infty)\\
and \\(-\infty,\infty)\\.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GPDDistrib.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GPDDistrib.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GPDDistrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GPDDistrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GPDDistrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GPDDistrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GPDDistrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.elementary.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GPDDistrib.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.GPDDistrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.GPDDistrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.GPDDistrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GPDDistrib.md).

Every other generic is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## Support depending on the parameters

The derivatives are correct as derivatives of the log-density at every
admissible point, whatever the sign of \\\xi\\. The license to
differentiate under the integral sign, on which the Bartlett identities
rest, is not automatic. The expected information exists for \\\xi \>
-1/2\\ and is Smith's closed form. At or below that value it does not
exist, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
returns `NA`; the classical asymptotics of the maximum likelihood
estimator do not hold either.

## See also

[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
to build one;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the \\\xi = 0\\ case;
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
for the companion family in an analysis of extremes;
[`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md)
for the upper limit of the support.

## Examples

``` r
d <- gpd_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

d@params
#> [1] "sigma" "xi"   
d@params_interpretation
#>   sigma      xi 
#> "scale" "shape" 

# The declared bounds are fixed at construction; the true endpoint moves
# with the parameters and is finite whenever the shape is negative.
c(declared = d@bounds,
  actual = c(0, distributions7:::gpd_endpoint(2, -0.4)))
#> declared1 declared2   actual1   actual2 
#>         0       Inf         0         5 

# The scale rides a log link and the shape may take either sign.
vapply(d@link_params, function(l) l@link_name, character(1))
#>      sigma         xi 
#>      "log" "identity" 
```
