# Dirichlet Distribution Class, Mean Vector and Concentration

The S7 class of the Dirichlet family on the open simplex of \\p\\
coordinates, parametrized by a mean vector \\\mu\\ carried on a
`parameters7` simplex and a concentration \\\phi \> 0\\. It inherits
from `multivariate_distrib`, so the distribution function and the
quantile are refused rather than approximated; the eleven methods listed
below are registered on it in this file and two more in `mv_higher.R`.

The class carries an extra property beyond the parent's, `param`: the
`parameters7` simplex the mean lies on, whose free names become the
family's own parameter names prefixed by `mean_`. Build one with
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md),
which validates the simplex against `n_dim`, supplies the
concentration's link and fills the properties in. This page documents
the raw S7 constructor, which validates none of that.

## Usage

``` r
DirichletDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  n_dim = integer(0),
  param = NULL
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

- n_dim:

  The dimension \\p\\ of one observation. A single positive integer; the
  validator rejects anything else, and also rejects a `dimension`
  property other than `"multivariate"`.

- param:

  The `parameters7`
  [`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
  carrying the mean. Its `n_free` is \\p-1\\, one fewer than the number
  of coordinates, the simplex being a set of dimension \\p-1\\ in
  \\\mathbb{R}^p\\.

## Value

An S7 object of class `DirichletDistrib`, inheriting from
`multivariate_distrib` and from `distrib`. Beyond `param` its properties
are the parent's: `distrib_name`, `dimension`, `n_dim`, `bounds`,
`params`, `params_interpretation`, `n_params`, `params_bounds`,
`link_params` and `params_smooth`. For an object built by
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
`params` is `c("mean_<free names>", "phi")`, `n_params` is \\p\\, and
every mean coordinate carries an identity link, its free value being
unconstrained already.

## Methods

Registered on this class in this file:
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.DirichletDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.DirichletDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.DirichletDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.DirichletDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.DirichletDistrib.md),
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.DirichletDistrib.md),
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.DirichletDistrib.md),
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.DirichletDistrib.md),
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.DirichletDistrib.md),
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.DirichletDistrib.md)
and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.DirichletDistrib.md).

Two more are registered in `mv_higher.R`:
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.DirichletDistrib.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.DirichletDistrib.md),
both closed form.

Everything else is inherited from
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md),
which refuses the distribution function, the quantile and the response
derivatives.

## See also

[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
to build one;
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md),
which is a coordinate's marginal and the two-coordinate case seen on the
line;
[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
for the discrete family on the same simplex;
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
for the chart.

## Examples

``` r
d <- dirichlet_distrib(3)
S7::S7_inherits(d, multivariate_distrib)
#> [1] TRUE

# Three coordinates, so two free mean values plus the concentration.
d@n_dim
#> [1] 3
d@params
#> [1] "mean_alr1" "mean_alr2" "phi"      
d@params_interpretation
#>       mean_alr1       mean_alr2             phi 
#>          "mean"          "mean" "concentration" 

# The mean lives on a parameters7 simplex, which is the extra property.
d@param@param_name
#> [1] "simplex"
d@param@n_free
#> [1] 2
```
