# S7 Class for Distributions With Fixed Parameters (Continuous)

The S7 class of a CONTINUOUS distribution in which some parameters of
the wrapped distribution are held at known values. It is the only
wrapper in this package that REMOVES parameters: the law is the
parent's, evaluated at a `theta` short by however many were fixed. Build
one with
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
which picks this class or one of its two siblings from the parent.

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
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

- parent_distrib:

  The wrapped `continuous_distrib` object.

- fixed_params:

  A named list of the fixed values, one single finite number each,
  strictly inside the corresponding parameter's open domain.

## Value

An S7 object of class `FixedContinuousDistrib`, inheriting from
`continuous_distrib` and from `distrib`. It carries `parent_distrib` and
`fixed_params` beside the parent's properties. For an object built by
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
`params` is the parent's less the fixed names, `n_params` falls by as
many, `bounds` is the parent's unchanged, and `distrib_name` is the
parent's with the held values in brackets.

## Methods

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
registers 22 methods on this class:
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md),
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md).

Every one splices the held values back into `theta` and delegates to the
parent. The derivative methods then subset the parent's answer by the
names the free parameter set generates, which is safe because the free
set preserves the parent's order: the same combination of free
parameters produces the same component string under either enumeration.

## How every method works

The free parameters are the parent's minus the fixed ones, IN THE
PARENT'S ORDER. Every method splices the fixed values back into `theta`
at their positions and delegates to the parent, so the parent's closed
forms are used wherever they exist, and a derivative method then keeps
only the components in which every index is a free parameter. No method
of this class computes anything of its own, and none is faster or more
accurate than the parent's.

## Why the components can be subset by name

The free set preserves the parent's ORDER, so a combination of free
parameters produces the same name string under the parent's enumeration
as under the wrapper's. Subsetting by name therefore cannot pair a
component with the wrong indices. That is the mistake re-parsing a name
by splitting on the underscore commits, for a parameter whose own name
contains one.

## The methods carry no documentation pages

The continuous and discrete branches register theirs inside a loop over
the two classes, so there is no top-level assignment for roxygen to
attach a block to. The multivariate branch registers its own at the top
level and is left undocumented to match, every one of them being the
same delegation. Read the parent's page for what a method computes and
this page for what the wrapper does to it.

## See also

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
to build one,
[FixedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/FixedDiscreteDistrib.md)
and
[FixedMultivariateDistrib](https://statmodels7.github.io/distributions7/reference/FixedMultivariateDistrib.md)
for the two siblings, and
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md),
whose half-normal is a fixed folded gaussian.

## Examples

``` r
# A gaussian centered at zero: one free parameter instead of two.
d <- fixed(gaussian1_distrib(), mu = 0)
d@params
#> [1] "sigma"
d@fixed_params
#> $mu
#> [1] 0
#> 
d@distrib_name
#> [1] "fixed gaussian1 [mu=0]"

# The law is the parent's at the spliced theta.
all.equal(distrib_pdf(d, c(-1, 0, 1), list(sigma = 2)),
          dnorm(c(-1, 0, 1), 0, 2))
#> [1] TRUE

# A derivative keeps only the components among the free parameters.
set.seed(1)
y <- rnorm(20, 0, 2)
names(distrib_gradient(d, y, list(sigma = 2)))
#> [1] "sigma"
names(distrib_deriv3(d, y, list(sigma = 2)))
#> [1] "sigma_sigma_sigma"

# Fixing everything is legal and gives a fully known law.
d0 <- fixed(gaussian1_distrib(), mu = 0, sigma = 1)
c(n_params = d0@n_params, density = distrib_pdf(d0, 0, list()))
#>  n_params   density 
#> 0.0000000 0.3989423 
```
