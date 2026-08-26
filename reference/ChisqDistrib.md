# Chi-Squared Distribution Class

The S7 class of the chi-squared family on \\(0, \infty)\\, parametrized
by its mean \\\mu \> 0\\, which is the degrees of freedom. They are
treated as a continuous positive parameter, so the one-parameter family
is estimable; the variance is then \\2\mu\\. It inherits from
`continuous_distrib`, so it answers every generic of the `distrib`
contract; the eleven methods listed below are registered on it directly
and everything else comes from the parent.

Build one with
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md),
which supplies the link function and fills the properties in. This page
documents the raw S7 constructor, which takes the parent's properties
and validates none of the relationships between them.

## Usage

``` r
ChisqDistrib(
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

An S7 object of class `ChisqDistrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
they hold `"chisq"`, `"univariate"`, `c(0, Inf)`, `"mu"`, the
interpretation `c(mu = "mean")`, `1`, the domain \\(0, \infty)\\, and
the one link.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ChisqDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ChisqDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ChisqDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ChisqDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ChisqDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ChisqDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ChisqDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ChisqDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ChisqDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ChisqDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ChisqDistrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
to build one;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
which contains this family at \\\sigma^2 = 2\mu\\;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the case \\\mu = 2\\;
[`distrib_pdf.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ChisqDistrib.md)
and
[`distrib_gradient.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ChisqDistrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- chisq_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# One parameter, and it is the mean.
d@params
#> [1] "mu"
d@n_params
#> [1] 1
d@params_interpretation
#>     mu 
#> "mean" 

# The variance is tied to it at 2 mu, so the family has no free spread.
vapply(c(1, 4, 20), function(m) variance(d, list(mu = m)), numeric(1))
#> [1]  2  8 40
```
