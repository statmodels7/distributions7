# S7 Class for Zero-Adjusted Discrete (Hurdle) Distributions

The S7 class of the hurdle version of a discrete distribution, with mass
function \$\$P(Y = y) = \begin{cases} \pi & y = 0 \\ (1-\pi)\dfrac{f(y;
\theta)}{1 - f(0; \theta)} & y \> 0. \end{cases}\$\$ The mass at zero is
REPLACED by the free parameter \\\pi\\, carried by `za`, and the parent
is truncated away from zero. Unlike
[ZeroInflatedDistrib](https://statmodels7.github.io/distributions7/reference/ZeroInflatedDistrib.md),
the model can produce FEWER zeros than the parent as well as more.

## Usage

``` r
ZeroAdjustedDiscreteDistrib(
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

An S7 object of class `ZeroAdjustedDiscreteDistrib`, inheriting from
`discrete_distrib` and from `distrib`. It carries `parent_distrib`
beside the parent's properties. For an object built by
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
the parameters are the parent's followed by `za`, whose bound is \\(0,
1)\\ and whose interpretation is `"prob. of zero"`.

## Details

The likelihood factorizes into a binary part in \\\pi\\ and a
positive-count part in \\\theta\\. Every mixed block of the Hessian is
exactly zero for that reason, and the two halves can be read separately.
It also re-interprets \\\theta\\: they are now the parameters of a
TRUNCATED law, not of the count process the observations come from.

Build one with
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
which checks that the parent is not already a zero wrapper and has
enough support points. This page documents the raw S7 constructor, which
checks neither.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroAdjustedDiscreteDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroAdjustedDiscreteDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroAdjustedDiscreteDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedDiscreteDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedDiscreteDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroAdjustedDiscreteDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ZeroAdjustedDiscreteDistrib.md)

The third and fourth derivatives come from the shared wrapper machinery
in `wrapper_derivatives.R`; everything else is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## Notation

\\f\\ is the parent's mass function, \\\pi\\ the probability of a zero
and \\\ell\\ the log-mass of one observation. The positive part is the
parent truncated away from zero, with mass \\f(y)/\\1 - f(0)\\\\ at \\y
\> 0\\.

## See also

[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
to build one,
[ZeroInflatedDistrib](https://statmodels7.github.io/distributions7/reference/ZeroInflatedDistrib.md)
for the wrapper that ADDS to the mass at zero, and
[ZeroAdjustedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedContinuousDistrib.md)
for the continuous branch, which produces a mixed distribution: a
density with an atom on it.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE
d@params
#> [1] "mu" "za"

# The mass at zero IS the parameter, and can be set below the parent's,
# which zero-inflation cannot do.
c(adjusted = distrib_pdf(d, 0, theta), parent = dpois(0, 3))
#>   adjusted     parent 
#> 0.40000000 0.04978707 
distrib_pdf(d, 0, list(mu = 3, za = 0.01))
#> [1] 0.01

# The positive part is the parent renormalized away from zero.
all.equal(distrib_pdf(d, 1:4, theta),
          0.6 * dpois(1:4, 3) / (1 - dpois(0, 3)))
#> [1] TRUE
```
