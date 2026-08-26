# S7 Class for Folded Distributions

The S7 class of the distribution of \\\|Y\|\\ when \\Y\\ follows the
wrapped continuous distribution, with density \$\$L(x; \theta) = f(x;
\theta) + f(-x; \theta), \qquad x \ge 0,\$\$ the two preimages of \\x\\
added together. It inherits from `continuous_distrib` and carries the
SAME parameters as its parent: folding adds none and removes none.

Build one with
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md),
which checks that the parent reaches below zero and carries no atom.
This page documents the raw S7 constructor, which validates neither.

## Usage

``` r
FoldedDistrib(
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

An S7 object of class `FoldedDistrib`, inheriting from
`continuous_distrib` and from `distrib`. It carries `parent_distrib`
beside the parent's `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
the parameters are the wrapped family's, `bounds` becomes
`c(0, max(abs(parent bounds)))`, and `distrib_name` is `"folded "`
followed by the parent's.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.FoldedDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.FoldedDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.FoldedDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.FoldedDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.FoldedDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.FoldedDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.FoldedDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.FoldedDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.FoldedDistrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md),
whose numerical quantile inverts the exact folded distribution function
above.

## See also

[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
to build one,
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the other wrapper that adds no parameter, and
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
which with `mu = 0` turns a folded gaussian into the half-normal.

## Examples

``` r
d <- folded(gaussian1_distrib())
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The parameters are the parent's, unchanged, and the support is the
# non-negative half line.
d@params
#> [1] "mu"    "sigma"
d@bounds
#> [1]   0 Inf
d@distrib_name
#> [1] "folded gaussian1"

# The parent is kept, so a method can reach it at both preimages.
d@parent_distrib@distrib_name
#> [1] "gaussian1"
```
