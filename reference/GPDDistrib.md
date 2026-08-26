# Generalized Pareto Distribution Class

The S7 class of the generalized Pareto family, the law of exceedances
over a high threshold, parametrized by a scale \\\sigma \> 0\\ and a
shape \\\xi\\: \\f(y) = \sigma^{-1}(1 + \xi y/\sigma)^{-1/\xi-1}\\ on
\\y \ge 0\\. At \\\xi = 0\\ it is the exponential, reached by a series,
so a fit may move the shape through zero.

It is the first family in the package whose **support moves with its
parameters**: for \\\xi \< 0\\ it ends at \\-\sigma/\xi\\. That is what
makes it non-regular, and it is why its expected information exists only
above \\\xi = -1/2\\.

Build one with
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md),
which supplies the two link functions. This page documents the raw S7
constructor, which validates none of the relationships between its
properties.

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

Registered in this file, all compiled apart from the first four:
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GPDDistrib.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GPDDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GPDDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GPDDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GPDDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GPDDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GPDDistrib.md).

Registered elsewhere: the third and fourth orders in `gpd_higher.R`
([`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GPDDistrib.md));
the response derivatives and the mixed one in
`cross_derivatives_families.R`
([`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md));
the four moments in `moments.R`; and the second-order response
derivatives in `theta2_more.R`.

## What the moving endpoint costs

The derivatives are correct as derivatives of the log-density at every
admissible point, whatever the sign of \\\xi\\. What is not automatic is
the license to differentiate under the integral sign, which the Bartlett
identities rest on. The expected information exists for \\\xi \> -1/2\\
and is Smith's closed form; at or below that it does not exist and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
returns `NA`, along with the classical asymptotics of the maximum
likelihood estimator.

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
