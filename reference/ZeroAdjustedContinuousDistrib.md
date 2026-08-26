# S7 Class for Zero-Adjusted Continuous Distributions

The S7 class of a continuous distribution with a point mass at zero:
\$\$P(Y = 0) = \pi, \qquad f_Y(y) = (1-\pi) f(y; \theta) \\ (y \ne
0).\$\$ The result is a MIXED distribution, a density plus an atom, and
it declares that atom through
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md).
That declaration is how
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
and
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
learn to treat it as one.

## Usage

``` r
ZeroAdjustedContinuousDistrib(
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

  The wrapped `continuous_distrib` object.

## Value

An S7 object of class `ZeroAdjustedContinuousDistrib`, inheriting from
`continuous_distrib` and from `distrib`. It carries `parent_distrib`
beside the parent's properties, with `za` added last.

## Details

No truncation is needed here. A continuous parent has \\P(Y = 0) = 0\\,
so there is no mass to remove before placing the atom, and the density
is simply scaled by \\1-\pi\\. That is the whole difference from
[ZeroAdjustedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedDiscreteDistrib.md),
whose parent must be truncated away from the point it already occupies.

The likelihood factorizes completely: the mixed blocks of the Hessian
are exactly zero, \\\pi\\ is estimated by the proportion of zeros and
\\\theta\\ by the parent's own fit to the non-zero observations.

Build one with
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).
This page documents the raw S7 constructor, which validates nothing.

## Methods

Registered on this class:
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroAdjustedContinuousDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroAdjustedContinuousDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ZeroAdjustedContinuousDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroAdjustedContinuousDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ZeroAdjustedContinuousDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedContinuousDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedContinuousDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroAdjustedContinuousDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ZeroAdjustedContinuousDistrib.md),
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.ZeroAdjustedContinuousDistrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## Notation

\\f\\ is the parent's density, \\\pi\\ the probability of a zero and
\\\ell\\ the log-density of one observation.

## See also

[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
to build one,
[ZeroAdjustedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedDiscreteDistrib.md)
for the discrete branch,
[`distrib_atoms.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md)
for the atom it declares, and
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md),
which REJECTS a parent of this class for exactly that atom.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)
d@params
#> [1] "mu"    "sigma" "za"   

# It is a MIXED distribution: an atom at zero and a density elsewhere.
distrib_atoms(d, theta)
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3
#> 
c(at_zero = distrib_pdf(d, 0, theta),
  elsewhere = distrib_pdf(d, 2, theta),
  scaled_parent = 0.7 * dnorm(2, 1, 2))
#>       at_zero     elsewhere scaled_parent 
#>     0.3000000     0.1232229     0.1232229 

# The density part integrates to 1 - pi, the atom carrying the rest.
integrate(function(z) ifelse(z == 0, 0, distrib_pdf(d, z, theta)),
          -Inf, Inf)$value
#> [1] 0.7
```
