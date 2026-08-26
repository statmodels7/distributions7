# Distribution Generics

The S7 generics every `distrib` object answers, grouped by what they
compute. This page is the index; each generic is documented on its own
page, and each family's method for it on a third.

**Only
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
is compulsory.** Everything else has a fallback: the distribution
function by quadrature or by summation, the quantile by root-finding or
by a table lookup, random draws by the ratio-of-uniforms method, and
every derivative by finite differences. A family defined by its density
alone therefore arrives complete, and a family that writes a quantity
out overrides the fallback by ordinary S7 dispatch.

## Value

Nothing. This page is the index of the generics documented on their own
pages; the value returned is theirs.

## The surface

- probability:

  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
  [`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md).

- derivatives in the parameters:

  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
  and
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md).

- derivatives in the response:

  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md),
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md).
  Refused for a discrete family: a lattice has nothing to differentiate
  along.

- mixed:

  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md),
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md).

- of the distribution function:

  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
  [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md),
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md),
  which a censored likelihood needs.

- moments:

  [`mean()`](https://rdrr.io/r/base/mean.html),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
  [`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md),
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md),
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md),
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

- fitting:

  [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md),
  and
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
  which is a function rather than a generic.

- multivariate:

  [`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md),
  [`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md),
  [`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
  [`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
  [`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md),
  [`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md).

## Two arguments every derivative generic takes

`theta` is normalized before dispatch by
[`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md),
which reorders it by name, strips stray names off the values and
validates against `params_bounds` treated as **open** intervals. `scale`
is applied by the generic **after** dispatch, so a method always returns
the parameter scale and never reads it; see
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

## See also

[`distrib()`](https://statmodels7.github.io/distributions7/reference/distrib.md)
for the base class and its properties;
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)
for the `scale` argument;
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the `approx` argument;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
which exercises the whole surface; the `defining-a-distribution`
vignette.

## Examples

``` r
# A family arrives complete from its density alone. Everything below the
# first line is a fallback the package supplies.
d <- gaussian1_distrib()
th <- list(mu = 0, sigma = 1)

distrib_pdf(d, 0, th)
#> [1] 0.3989423
distrib_cdf(d, 0, th)
#> [1] 0.5
distrib_quantile(d, 0.975, th)
#> [1] 1.959964
distrib_gradient(d, c(-1, 1), th)
#> $mu
#> [1] -1  1
#> 
#> $sigma
#> [1] 0 0
#> 

# theta is reordered by name before dispatch, so the order it is written in
# does not matter.
identical(distrib_pdf(d, 1, list(mu = 0, sigma = 2)),
          distrib_pdf(d, 1, list(sigma = 2, mu = 0)))
#> [1] TRUE

# The bounds are OPEN, so a scale of exactly zero is rejected rather than
# returning an infinite density.
try(distrib_pdf(d, 1, list(mu = 0, sigma = 0)))
#> Error : Invalid parameter value(s) for the 'gaussian1' distribution:
#>   'sigma' = 0 is outside its domain (0, Inf)

# A discrete family refuses the response derivatives by design.
try(distrib_grad_y(poisson_distrib(), 2, list(mu = 3)))
#> Error : Can't find method for `distrib_grad_y(<distributions7::PoissonDistrib>)`.
```
