# Mixed Response Derivatives of a Zero-Adjusted Continuous Family

Evaluates a mixed derivative in the response and the parameters,
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
or
[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md),
for
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
over a continuous parent. Away from the atom the components of the
parent's parameters are the parent's, and the component of `za` is zero,
since `za` enters the log-density only through \\\log(1 - \pi)\\. At the
atom every component is `NaN`, as for
[`za_y_deriv()`](https://statmodels7.github.io/distributions7/reference/za_y_deriv.md).
The same guard is registered for
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
and
[`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
through
[`za_y_deriv()`](https://statmodels7.github.io/distributions7/reference/za_y_deriv.md).

## Usage

``` r
za_y_list(distrib, y, theta, fun)
```

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with the parent's parameters followed by `za`.

- fun:

  The generic to evaluate on the parent.

## Value

A named list keyed by `distrib@params`, one numeric vector per
parameter, each as long as `y`.

## See also

[`za_y_deriv()`](https://statmodels7.github.io/distributions7/reference/za_y_deriv.md)
for the response derivatives.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
distributions7:::za_y_list(d, c(-1, 0, 2), list(mu = 1, sigma = 2, za = 0.3),
                           distrib_cross_y)
#> $mu
#> [1] 0.25  NaN 0.25
#> 
#> $sigma
#> [1] -0.50   NaN  0.25
#> 
#> $za
#> [1]   0 NaN   0
#> 
```
