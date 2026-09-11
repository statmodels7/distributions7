# Numerically Validate a Distribution

Runs a battery of numerical self-consistency checks on a `distrib`
object. Validates a user-defined distribution: it verifies that the
density integrates (or sums) to one, that the CDF, quantile function and
random generator agree with each other and with the density, and that
every analytical derivative matches its finite-difference counterpart.

## Usage

``` r
check_distrib(
  distrib,
  theta = NULL,
  n = 100,
  nsim = 2e+05,
  orders = 1:4,
  tol = 0.001,
  verbose = TRUE
)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameter values at which to run the checks. If `NULL`
  (default) a random admissible value is drawn with
  [`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md).

- n:

  Integer. Number of observations used for the derivative comparisons.
  Defaults to 100.

- nsim:

  Integer. Monte Carlo sample size used for the random generator and
  expected-information checks. Defaults to 200000.

- orders:

  Integer vector. Which parameter-derivative orders to check. Defaults
  to `1:4`, the orders every family implements analytically; use e.g.
  `1:2` for a faster run. Add `5` to check the numerical fifth as well;
  see the bullet below for what that row compares and when it is
  emitted.

- tol:

  Numeric. Relative tolerance for the finite-difference comparisons.
  Defaults to `1e-3`.

- verbose:

  Logical. If `TRUE` (default) a readable report is printed.

## Value

Invisibly, a `data.frame` with one row per check and columns `check`,
`status` (`"OK"` or `"FAIL"`), `statistic` and `detail`. A check whose
statistic came out `NaN` or `NA` is not a row: it is listed in the
attribute `"skipped"`, a data frame with columns `check` and `reason`,
which is present only when some check was skipped.

## Details

The checks performed are:

- **density**: non-negativity and integration/summation to 1 over the
  support.

- **cdf**: values in \\\[0,1\]\\ and monotonicity along a grid of
  quantiles, and its agreement with the density: for a continuous family
  a central difference of the cdf against the density, with the step the
  response derivatives use, and for a discrete one \\F(k) - F(k-1)\\
  against the mass.

- **quantile**: round-trip against the CDF (\\F(Q(p)) = p\\ for
  continuous distributions, and the generalized-inverse inequalities for
  discrete ones).

- **rng**: the sample mean and variance of a large draw agree with
  [`mean()`](https://rdrr.io/r/base/mean.html) and
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
  within Monte Carlo error.

- **gradient, hessian, deriv3, deriv4**: analytical values against
  [`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md),
  [`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md),
  [`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md)
  and
  [`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md).

- **deriv5**, when `5` is among `orders`:
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
  against
  [`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md)
  at a higher accuracy, five stencil nodes instead of three. No family
  writes the fifth order out, so there is no analytic value to compare
  against and this checks the differencing rather than a family's
  algebra. It is emitted only where
  [`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
  is `TRUE`: where the fourth order is itself a fallback the fifth is a
  difference of a difference and no verdict on it would mean anything. A
  family that owns its fourth-order method while building part of it
  from stencils passes that test and may still fail this row:
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
  fails it at `nu = 3` and passes from `nu = 8` upward, the noise read
  here being absolute and falling as `nu` grows. It is why `orders`
  defaults to `1:4`.

- **expected information**:
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
  against a Monte Carlo estimate of
  \\-\mathbb{E}\[\nabla\ell\\\nabla\ell^\top\]\\. The outer product of
  the score is used as reference because it remains valid when the
  log-likelihood is not differentiable in a parameter (see
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)).
  A draw landing exactly on a finite bound of a continuous support,
  where the score of a family singular at that bound is infinite, is
  left out of the estimate and counted in the row's `detail`; a score
  that is not finite at an interior point still leaves the row without a
  statistic.

- **response derivatives** (continuous only):
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
  and
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
  against central differences in \\y\\ whose step is chosen observation
  by observation by
  [`fd_stable_quotient()`](https://statmodels7.github.io/distributions7/reference/fd_stable_quotient.md),
  between a step cut to under half the distance to a bound and one
  scaled on that distance, each divided by the steps its evaluation
  points actually lie at. A draw closer than \\128\lvert
  b\rvert\varepsilon\\ to a finite bound \\b \ne 0\\, where the spacing
  of doubles is absolute and no central difference compares a
  derivative, is left out and counted in the row's `detail`.

- **link scale**: `scale = "link"` derivatives against finite
  differences of the log-likelihood in \\\eta\\.

Distributions that rely on the numerical fallbacks pass the
corresponding parameter-derivative checks trivially, since analytical
and numerical values then coincide by construction. The response
fallbacks take the cut step alone, so near a bound they can differ from
the reference, which is what the row then reports.

A check whose statistic comes out `NaN` or `NA` has nothing to judge, as
the expected information of a family where it does not exist, and is not
a row of the table: it is listed with its reason in the attribute
`"skipped"`, and every row that is in the table reads `"OK"` or
`"FAIL"`. An infinite statistic is a failure, being a component that
overflows where its reference does not; so is a check whose computation
raised an error, and so is a density or a distribution function that is
not finite, those three checks being defined on the values themselves.

Mixed distributions — a density with point masses on top of it, as
produced by
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
on a continuous parent — are handled as long as they declare their atoms
through
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md).
The density is then expected to integrate to one minus the atomic mass,
quantiles falling inside a jump of the CDF are checked as generalized
inverses rather than exact ones, and finite differences in \\y\\ are
kept away from the atoms, where no derivative exists.

## See also

[`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md),
[`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md),
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)

## Examples

``` r
if (FALSE) { # \dontrun{
check_distrib(gaussian1_distrib())
check_distrib(laplace_distrib(), theta = list(mu = 1, sigma = 2))
check_distrib(poisson_distrib(), orders = 1:2, nsim = 5e4)
} # }
```
