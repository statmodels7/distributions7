# S7 Class for Truncated Continuous Distributions

Represents a continuous parent restricted to \\\[L, U\]\\ and
renormalized by the retained mass \\Z(\theta)\\. Construct one with
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
which validates the endpoints, collapses a nested truncation and copies
the parent's parameter metadata; calling the class directly does none of
that.

## Usage

``` r
TruncatedContinuousDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  parent_distrib = distrib(),
  lower = integer(0),
  upper = integer(0)
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

  The wrapped `continuous_distrib` object. Its parameters become the
  truncated object's, unchanged.

- lower, upper:

  The truncation points, `L` and `U`. Either may be infinite, giving
  one-sided truncation, and both are included in the support.

## Value

An S7 object of class `TruncatedContinuousDistrib`, inheriting from
`continuous_distrib`.

## What truncation adds, and what it does not

It adds NO parameter. The endpoints are known constants, like a
binomial's `size`, so the truncated object carries exactly the parent's
`params`, `params_bounds` and `link_params`. What it adds is the
\\\theta\\-dependent normalizing constant \\Z\\, and every derivative of
\\\ell_T = \ell - \log Z\\ carries its contribution.

The support does not move with \\\theta\\. That is the condition under
which \\Z\\ may be differentiated under the integral sign, and it is
what keeps truncation at fixed points a regular problem.

## A mixed parent

The parent may itself carry point masses, as
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
of a continuous distribution does. Those atoms survive truncation where
they lie inside the interval, rescaled by \\1/Z\\, and
[`distrib_atoms.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.TruncatedContinuousDistrib.md)
reports them. Two methods exist only for that case:
[`expectation.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/expectation.TruncatedContinuousDistrib.md)
adds the masses to the integral, and
[`parent_mass_at()`](https://statmodels7.github.io/distributions7/reference/parent_mass_at.md)
asks the parent for the mass on a single point instead of assuming a
continuous parent has none.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## Methods

Methods implemented for this class:
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.TruncatedContinuousDistrib.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedContinuousDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedContinuousDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.TruncatedContinuousDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TruncatedContinuousDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.TruncatedContinuousDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedContinuousDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedContinuousDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedContinuousDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedContinuousDistrib.md),
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.TruncatedContinuousDistrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the constructor,
[TruncatedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedDiscreteDistrib.md)
for the discrete branch, and
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
for \\Z\\.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

class(tn)[1]
#> [1] "distributions7::TruncatedContinuousDistrib"
c(name = tn@distrib_name, params = paste(tn@params, collapse = ", "))
#>                          name                        params 
#> "truncated gaussian1 [-1, 2]"                   "mu, sigma" 
tn@bounds
#> [1] -1  2

# The parameters are the parent's, unchanged: truncation adds none.
identical(tn@params, gaussian1_distrib()@params)
#> [1] TRUE
identical(tn@params_bounds, gaussian1_distrib()@params_bounds)
#> [1] TRUE

# The density is the parent's divided by the retained mass.
Z <- pnorm(2, 0.3, 1.2) - pnorm(-1, 0.3, 1.2)
c(truncated = distrib_pdf(tn, 0, theta),
  parent_over_Z = dnorm(0, 0.3, 1.2) / Z)
#>     truncated parent_over_Z 
#>     0.4118505     0.4118505 
```
