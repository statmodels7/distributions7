# S7 Class for Zero-Inflated Distributions

The S7 class of the zero-inflated version of a discrete distribution,
with mass function \$\$P(Y = y) = \zeta\\\mathbb{I}(y = 0) + (1 -
\zeta)\\ f(y; \theta).\$\$ It inherits from `discrete_distrib` and
carries the parent's parameters followed by `zi`, which is the
probability \\\zeta\\ of a structural zero and rides a link of its own.

Build one with
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
which checks that the parent is discrete, is not already a zero wrapper,
and has enough support points for the extra parameter to be identified.
This page documents the raw S7 constructor, which checks none of that.

## Usage

``` r
ZeroInflatedDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  parent_distrib = distrib()
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

- parent_distrib:

  The wrapped `discrete_distrib` object.

## Value

An S7 object of class `ZeroInflatedDistrib`, inheriting from
`discrete_distrib` and from `distrib`. It carries `parent_distrib`
beside the parent's `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
the parameters are the parent's followed by `zi`, whose bound is \\(0,
1)\\ and whose interpretation is `"prob. of structural zero"`.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroInflatedDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroInflatedDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroInflatedDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroInflatedDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroInflatedDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroInflatedDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ZeroInflatedDistrib.md)

The third and fourth derivatives come from the shared wrapper machinery
in `wrapper_derivatives.R`; everything else is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## See also

[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
to build one,
[ZeroAdjustedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedDiscreteDistrib.md)
for the wrapper that REPLACES the mass at zero where this one adds to
it, and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
the constructor to reach for with a continuous parent.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

# The parent's parameters, then zi last with its own link and bound.
d@params
#> [1] "mu" "zi"
d@params_bounds$zi
#> [1] 0 1
vapply(d@link_params, function(l) l@link_name, character(1))
#>      mu      zi 
#>   "log" "logit" 
d@params_interpretation
#>                         mu                         zi 
#>                     "mean" "prob. of structural zero" 

# Inflation can only ADD zeros: the mass at zero exceeds the parent's.
c(inflated = distrib_pdf(d, 0, theta), parent = dpois(0, 3))
#>   inflated     parent 
#> 0.28734030 0.04978707 
```
