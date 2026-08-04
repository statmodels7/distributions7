# S7 Class for Distributions With Fixed Parameters (Continuous)

A subclass of `continuous_distrib` representing a continuous
distribution in which some parameters of the wrapped distribution are
held at known values. Constructed by
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.md).

## Usage

``` r
FixedContinuousDistrib(
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
  fixed_params = list()
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
  [`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

- parent_distrib:

  The wrapped `continuous_distrib` object.

- fixed_params:

  A named list of the fixed parameter values.

## Value

An object of class `FixedContinuousDistrib`.

## Details

The free parameters are the parent's minus the fixed ones, in the
parent's order. Every method splices the fixed values back into `theta`
at their positions and delegates to the parent, so the parent's closed
forms are used whenever they exist; a derivative method then keeps only
the components in which every index is a free parameter. No method of
this class computes anything of its own.

## See also

[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.md)
