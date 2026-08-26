# Print Method for `distrib` Objects

Writes what a distribution object is: its name and dimension, the
support, and one row per parameter giving the parameter's name, what it
means, the interval it lives in and the link that carries it to the
unconstrained scale. A wrapper prints its parent's line first, so
`truncated(gamma2_distrib(), upper = 5)` shows both the truncation and
the family under it.

What it does **not** show is a value: a `distrib` carries a
parametrization, not an estimate. For estimates see
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md).

## Arguments

- x:

  An object inheriting from `distrib`.

- ...:

  Unused, accepted for compatibility with
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly. Called for the output it writes.

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
