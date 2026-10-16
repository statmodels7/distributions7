# Print Method for `distrib` Objects

Prints the name of the distribution, its type (continuous or discrete)
and its dimension, followed by one row per parameter giving the name of
the parameter, its interpretation, the interval in which it lies and the
link that maps it to the unconstrained scale. The name of a wrapper
combines the wrapper and its parent, so
`truncated(gamma2_distrib(), upper = 5)` prints as
`Truncated Gamma2 [upper=5]`, followed by the parameter rows of the
parent. No parameter value is shown, because a `distrib` carries a
parametrization and no estimate; for estimates see
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md).

## Arguments

- x:

  An object inheriting from `distrib`.

- ...:

  Unused, accepted for compatibility with
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly. The function is called for the output that it writes.

## See also

[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
for a fitted object;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate one;
[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
to draw one.

## Examples

``` r
gaussian1_distrib()
#> Distribution: Gaussian1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   sigma (standard deviation) | Link: log        | Domain: (0, Inf)
poisson_distrib()
#> Distribution: Poisson
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (mean)               | Link: log        | Domain: (0, Inf)

# The link is part of the object, so changing it changes what prints.
poisson_distrib(link_mu = linkfunctions7::sqrt_link())
#> Distribution: Poisson
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (mean)               | Link: sqrt       | Domain: (0, Inf)
```
