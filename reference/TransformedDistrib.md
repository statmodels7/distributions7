# S7 Class for Transformed Distributions

The S7 class of the distribution of \\Y = g(X)\\, where \\X\\ follows a
wrapped continuous distribution and \\g\\ is a BIJECTIVE
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md).
Its density is the change-of-variables formula \$\$f_Y(y) =
f_X(g^{-1}(y))\\\lvert J(y)\rvert, \qquad J(y) =
\frac{d}{dy}g^{-1}(y),\$\$ and it carries exactly the parent's
parameters: the transformation adds none and removes none.

## Usage

``` r
TransformedDistrib(
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
  transformer = transformer()
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

- transformer:

  The
  [`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
  defining \\g\\.

## Value

An S7 object of class `TransformedDistrib`, inheriting from
`continuous_distrib` and from `distrib`. It carries `parent_distrib` and
`transformer` beside the parent's properties. For an object built by
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
`params`, `params_bounds` and `link_params` are the parent's unchanged;
`bounds` is the transformer's image of the parent's;
`params_interpretation` is the parent's with the parent's name appended
in brackets, since a mean on the transformed scale is not the parent's
mean; and `distrib_name` is `"g(parent)"` or whatever `new_name` said.

## Why the derivatives are the parent's

\\g\\ does not depend on \\\theta\\, so the Jacobian is a constant in
the parameters and \\\ell_Y(\theta; y) = \ell_X(\theta; g^{-1}(y)) +
\log\lvert J(y)\rvert\\ differentiates in \\\theta\\ to the parent's own
derivative at \\x = g^{-1}(y)\\. The score, the observed Hessian and the
EXPECTED Hessian are therefore the parent's, exactly, with no Monte
Carlo anywhere. The response derivatives are not: those pick the
Jacobian up and come from the base class.

Build one with
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which checks that the parent is continuous and that the transformer is
valid on its support. This page documents the raw S7 constructor, which
checks neither.

## Notation

\\g\\ is the transformation, \\J\\ the Jacobian of its inverse, \\f_X\\
the parent's density, \\f_Y\\ the transformed one and \\\ell\\ a
log-density.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TransformedDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TransformedDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TransformedDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TransformedDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TransformedDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TransformedDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TransformedDistrib.md)

The third and fourth derivatives come from the shared wrapper machinery
in `wrapper_derivatives.R`; everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
to build one,
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
for the map,
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
for a map that is NOT bijective and so cannot be a transformer, and
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
for the wrapper that removes parameters.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The parameters are the parent's, and the support is the image of its own.
d@params
#> [1] "mu"    "sigma"
d@bounds
#> [1]   0 Inf
d@distrib_name
#> [1] "exp(gaussian1)"

# The exponential of a gaussian is the lognormal, whose second parameter is
# the variance rather than the standard deviation.
y <- c(0.5, 1, 3)
all.equal(distrib_pdf(d, y, theta),
          distrib_pdf(lognormal1_distrib(), y,
                      list(mu = 0.5, sigma2 = 0.8^2)))
#> [1] TRUE

# The interpretation says which scale a parameter lives on.
d@params_interpretation
#>                                     mu                                  sigma 
#>               "mean (gaussian1 scale)" "standard deviation (gaussian1 scale)" 
```
