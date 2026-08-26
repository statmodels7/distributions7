# S7 Class for Truncated Discrete Distributions

Represents a discrete parent restricted to the support points in \\\[L,
U\]\\ and renormalized by the retained mass \\Z(\theta)\\. Both
endpoints are INCLUDED, so `truncated(poisson_distrib(), lower = 1)` is
the zero-truncated Poisson on \\\\1, 2, \dots\\\\ and keeps the mass at
one. Construct one with
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md);
calling the class directly skips every validation the constructor
performs.

## Usage

``` r
TruncatedDiscreteDistrib(
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

  The wrapped `discrete_distrib` object. Its parameters become the
  truncated object's, unchanged.

- lower, upper:

  The truncation points, `L` and `U`, both included in the support.
  Either may be infinite, and a finite one must be a whole number.

## Value

An S7 object of class `TruncatedDiscreteDistrib`, inheriting from
`discrete_distrib`.

## Details

The lower endpoint being included is the one place the two truncation
classes differ. The tail below the interval is \\F(L^-) = F(L) - f(L)\\,
so the mass sitting exactly on \\L\\ has to be added back;
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
does that through
[`parent_mass_at()`](https://statmodels7.github.io/distributions7/reference/parent_mass_at.md),
and
[`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md)
applies the matching correction to the derivatives of \\Z\\.

Truncation adds no parameter, so the object carries exactly the parent's
`params`, `params_bounds` and `link_params`. What it can remove is
IDENTIFIABILITY: \\k\\ retained support points carry \\k-1\\ free
probabilities, so
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
rejects an interval leaving fewer than `n_params + 1` of them.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## Methods

Methods implemented for this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedDiscreteDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedDiscreteDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TruncatedDiscreteDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedDiscreteDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedDiscreteDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedDiscreteDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedDiscreteDistrib.md)

Everything else is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the constructor,
[TruncatedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedContinuousDistrib.md)
for the continuous branch, and
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
for \\Z\\.

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
class(ztp)[1]
#> [1] "distributions7::TruncatedDiscreteDistrib"
c(name = ztp@distrib_name, bounds = paste(ztp@bounds, collapse = ", "))
#>                          name                        bounds 
#> "truncated poisson [lower=1]"                      "1, Inf" 

# The lower endpoint is IN the support, so this is the zero-truncated
# Poisson: mass at 0 is gone, mass at 1 is not.
distrib_pdf(ztp, 0:3, list(mu = 2))
#> [1] 0.0000000 0.3130353 0.3130353 0.2086902
dpois(1:3, 2) / (1 - dpois(0, 2))
#> [1] 0.3130353 0.3130353 0.2086902

# An interval leaving too few support points is rejected by the
# constructor, the parameter no longer being identified.
try(truncated(bernoulli_distrib(), lower = 1))
#> Error : lower = 1 lies at or above the support of 'bernoulli', which ends at 1: the truncated support would be empty.
```
