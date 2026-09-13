# Numerical Mixed Third-Response Parameter Derivatives

Computes \\\partial^4 \ell / \partial y^3\\ \partial \theta_i\\ by one
central difference of
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
in each parameter. It is what the default
[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md)
method runs for a continuous family that is not a location family.

## Usage

``` r
numerical_cross3_y(
  distrib,
  y,
  theta,
  h_rel = .Machine$double.eps^(1/3),
  which = NULL
)
```

## Arguments

- distrib:

  A distribution object inheriting from `continuous_distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, the point to differentiate at.

- h_rel:

  The relative step, defaulting to `.Machine$double.eps^(1/3)`, the
  optimal exponent for a central first difference.

- which:

  An optional character vector naming a subset of `distrib@params`. The
  default, `NULL`, computes every component.

## Value

A named list with one numeric vector per requested parameter, each as
long as `y`.

## Details

A family carrying an analytic `distrib_deriv3_y` pays for exactly one
difference. Where that derivative is itself the base class's stencil on
the log-density, the two differences act on different variables and
compose into one mixed stencil. The step is chosen by
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md),
so a parameter near a bound is differenced on the side the bound allows.

## See also

[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md),
the generic it serves, and
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
for the quantity it differences.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3, nu = 6)

# Against the identity the location family uses.
max(abs(unlist(numerical_cross3_y(d, y, theta)) -
        unlist(distrib_cross3_y(d, y, theta))))
#> [1] 8.227036e-11

numerical_cross3_y(d, y, theta, which = "nu")
#> $nu
#> [1]  0.03085952  0.02819720 -0.02265895
#> 
```
